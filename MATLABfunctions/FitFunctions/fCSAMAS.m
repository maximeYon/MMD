function Y = fCSAMAS(Pin,Xin,Pnorm,Xnorm,Ynorm,Para);

Pin = Pin.*Pnorm;
nu = Xin*Xnorm;

Intensity = Pin(1);
T2 = Pin(2);
deltaiso = Pin(3);
deltaaniso = Pin(4);
eta = Pin(5);

nu0 = Para(1);
nur = Para(2);

deltazz = deltaiso + deltaaniso;
deltayy = (3*deltaiso - eta*deltaaniso - deltazz)/2;
deltaxx = eta*deltaaniso + deltayy;

if deltaaniso < 0
    minSB = round(deltazz/(nur/nu0*1e6));
    maxSB = round(deltayy/(nur/nu0*1e6));
else
    maxSB = round(deltazz/(nur/nu0*1e6));
    minSB = round(deltayy/(nur/nu0*1e6));
end

RefAxis = (minSB-3):(maxSB+3);

qu = 6;
NrGammaAngles=max([50 2*max(abs(RefAxis))+10]);

%
%Convert the parameters to their proper units
%
omega0 = nu0*pi*2; 		%2*pi*nu with nu in Hz
omegar = nur*pi*2; 		%2*pi*nu with nu in Hz

nsteps = NrGammaAngles;
tresolution = 2*pi/omegar/nsteps;
tarray = 0:tresolution:(nsteps-1)*tresolution;

%
%perform the gamma-averaged calculation (follows the symvols from the gamma-COMPUTE paper!)
%

%
%values for number of points in two-angle powder averaging and the weighting function
%
pa = [3, 21, 144, 233, 399, 610, 987, 1583, 2741, 4409, 6997, 11657, 17389, 28499, 43051, 65063, 79999, 2943631];
pb = [2, 13,  37, 163, 359, 269, 937, 1153, 1117, 1171, 3049,   131,  1787,  5879,  4649,  5237, 74729,  702707];
NrPP = pa(qu);
beta  = pi * ((1:NrPP)-1)/NrPP;
alpha = 2*pi * (mod(pb(qu)*((1:NrPP)-1), NrPP)/NrPP);
gamma = 0;
weight = sin(beta) / NrPP;
weightMatrix = weight.' * ones(1, nsteps);


%
%determine the sideband manifold requested by the user through the RefAxis
%
FFTstart = floor((-nsteps+1)/2); 	%starting point of the sideband manifold
FFTstop  = nsteps -1 + FFTstart; 	%end point of the sideband manifold
FFTNrSidebands = nsteps;
FFTvector = FFTstart:FFTstop;

for tel=1:length(RefAxis)
  SSAindexFFT(tel) = find(FFTvector == RefAxis(tel));
end
SSA = 0*RefAxis;


%
%calculate all frequencies for all crystallites at all times.
%
omegars = FrequencyMAS(omega0, deltaaniso*1e-6, eta, omegar, alpha, beta, gamma, tarray);

%figure(1), clf, plot(tarray',omegars'), return

%
%these are the QTrs factors for all crystallites and all time steps
%
QTrs = ones(NrPP, nsteps);
for j=2:nsteps
  QTrs(:, j) = exp(-sqrt(-1)*sum(omegars(:, 1:(j-1)), 2)*tresolution);
end


%
%these are the rhoT0sr factors for all crystallites and all time steps
%
rhoT0sr = conj(QTrs);


%
%calculate the gamma-averaged FID over 1 rotor period for all crystallites
%
for j=0:nsteps-1
  favrs(:, j+1) = sum(rhoT0sr .* (QTrs(:, mod( (0:nsteps-1)+j, nsteps) +1)), 2);
end


%
%apply the weight of the powder averaging scheme (sideband intensities are independent of nsteps and NrPP!!)
%
favrs = favrs .* weightMatrix / nsteps^2;

%figure(1), clf, plot(tarray',favrs'), return
%figure(1), clf, plot(tarray',sum(favrs,1)'), return

%
%calculate the sideband intensities by doing an FT and pick the ones that are needed further
%
sidebands = fftshift(real(fft(sum(favrs))));
sidebands = sidebands/sum(sidebands);
%figure(1), clf, plot(sidebands), return
SSA = Intensity*sidebands(SSAindexFFT);

%figure(1), clf, plot(RefAxis,SSA,'o')

ampl = SSA;
offset = RefAxis*nur/nu0*1e6 + deltaiso;

I = fSumLorentz([T2*ones(size(SSA)) ampl offset],nu,1,1,1);

% figure(3), clf
% plot(nu,I,'-')
% hold on
% set(gca,'XDir','reverse')
% plot([deltaxx deltayy deltazz],[0 0 0],'x')

Y = I/Ynorm;
