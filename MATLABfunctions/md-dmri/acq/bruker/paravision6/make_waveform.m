% Gradient waveforms for axisymmetric b-tensors.
%
% Calculated according to the recipe in fig 1 of
% Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016).
% http://dx.doi.org/10.1039/c5cp07251d
%
% Output in the format for Matthew Budde's code for Bruker pv6.0.1 

% Define path for output folders
run_path = fileparts(mfilename('fullpath'));
out_path = fullfile(run_path,'waveforms');


% Parameters for calculating b-values
gmax = 100/100*0.7853; % Max gradient [3 T/m]
tau = 15e-3; % Waveform duration [10 ms]

% Define timing parameters relative to a total echo time = 1
np = 500; % Number of time intervals in waveform [1000]
epsilon_up = .1; % Quarter-sine ramp up [0.1]
plateau = 0; % Plateau [0]
epsilon_down = .15; % Half-sine ramp down [0.1]

% q-trajectory parameters
% zeta: half aperture of q cone, see Topgaard. Microporous Mesoporous Mater. 178, 60 (2013).
% http://dx.doi.org/10.1016/j.micromeso.2013.03.009
% zeta = acos(1/sqrt(3)); out_fn = fullfile(out_path,['axde_sphere_' num2str(np)]); %sphere
% zeta = pi/2; out_fn = fullfile(out_path,['axde_plane_' num2str(np)]); %plane
% zeta = 0; out_fn = fullfile(out_path,['axde_stick_' num2str(np)]); %stick
% zeta = acos(sqrt(2/3)); out_fn = fullfile(out_path,['axde_cigar_' num2str(np)]); %cigar
zeta = -1; out_fn = fullfile(out_path,['axde_stick_hf_' num2str(np)]); %stick high-frequency

% (theta,phi): Orientation of cone axis in lab frame
theta = 0; phi = 0;

%------------------------

taured = 1;
dt = taured/np;
t = taured*linspace(0,1,np);

% Quarter-sine ramp up
np_epsilon_up = round(epsilon_up/dt);
t_epsilon_up = pi/2*linspace(0,1,np_epsilon_up)';
g_up = sin(t_epsilon_up);
%figure(1), clf, plot(t_epsilon_up,g_up,'-'), return

% Half-sine ramp down
np_epsilon_down = round(epsilon_down/dt);
t_epsilon_down = pi/2*linspace(-1,1,np_epsilon_down)';
g_down = 1 - .5*(1+sin(t_epsilon_down));
%figure(1), clf, plot(t_epsilon_down,g_down,'-'), return

% Plateau
np_plateau = round(plateau/dt);
np_inter = round(taured/dt)-2*(np_epsilon_up+np_plateau+np_epsilon_down);

ga = [g_up; ones(np_plateau,1); g_down; zeros(np_inter,1); g_down-1; -1*ones(np_plateau,1); -flipud(g_up)];
%figure(1), clf, plot(t,ga,'-'), return

deltapsi = 2*pi;
psi0 = 0;
b_delta = (3*cos(zeta).^2-1)/2; % normalized anisotropy of the b tensor

q = cumsum(ga); % Topgaard Eq 35 (Note: all dt cancel and can be omitted.)
q = mean([q flipud(q)],2); % Removes asymmetry introduced by cumsum
b = sum(q.^2); % Topgaard Eq 32
psi = psi0 + deltapsi/b*cumsum(q.^2); % Topgaard Eq 31
psi = mean([psi deltapsi-flipud(psi)],2); % Removes asymmetry introduced by cumsum
gr = (ga + 1i*deltapsi/b*q.^3).*exp(1i.*psi); % Topgaard Eq 37

% Assure that first and last values equal zero
ga([1 end]) = 0;
gr([1 end]) = 0;

re_gr = real(gr);
im_gr = imag(gr);

% Assure refocusing on each channel
ga(2:(end-1)) = ga(2:(end-1)) - sum(ga)/(np-2);
re_gr(2:(end-1)) = re_gr(2:(end-1)) - sum(re_gr)/(np-2);
im_gr(2:(end-1)) = im_gr(2:(end-1)) - sum(im_gr)/(np-2);

% Normalize waveform to the range -1 to +1
gnorm = max(abs([ga(:); re_gr(:)+1i*im_gr(:)]));
ga = ga/gnorm;
re_gr = re_gr/gnorm;
im_gr = im_gr/gnorm;
gr = re_gr+1i*im_gr;

gx = re_gr*sin(zeta); % Topgaard Eq 37
gy = im_gr*sin(zeta);
gz = ga*cos(zeta);

if zeta == -1
    gx = 0*re_gr;
    gy = 0*im_gr;
    gz = re_gr;
end

gx_old = gx;
gz_old = gz;
gx = gx_old*cos(theta) + gz_old*sin(theta);
gz = gz_old*cos(theta) - gx_old*sin(theta);
gx_old = gx;
gy_old = gy;
gx = gx_old*cos(phi) - gy_old*sin(phi);
gy = gy_old*cos(phi) + gx_old*sin(phi);

gamma = 26.75e7;
dt = tau/np;
t = tau*linspace(0,1,np)';

qx = gamma*gmax*cumsum(gx*dt);
qy = gamma*gmax*cumsum(gy*dt);
qz = gamma*gmax*cumsum(gz*dt);
% qx = cumsum(gx);
% qy = cumsum(gy);
% qz = cumsum(gz);
bxx = sum(qx.^2);
byy = sum(qy.^2);
bzz = sum(qz.^2);
bxy = sum(qx.*qy);
bxz = sum(qx.*qz);
byz = sum(qy.*qz);

b_tensor_norm = [[bxx bxy bxz]; [bxy byy byz]; [bxz byz bzz]]/(bxx+byy+bzz)
% If isotropic, normalized b tensor is diagonal with all eigenvalues 1/3

% q = gamma*gmax*cumsum(ga*dt);
q = sqrt(qx.^2 + qy.^2 + qz.^2);
b = sum(q.^2*dt);

qmax = max(q);
td = b/qmax^2;
bfactor = b/gmax.^2;

dgadt = gradient(gmax*ga/dt);
dgrdt = gradient(gmax*abs(gr)/dt);
dre_grdt = gradient(gmax*re_gr/dt);
dim_grdt = gradient(gmax*im_gr/dt);

figure(1), clf
subplot(2,2,1)
%plot(t,gmax*re_gr,'r-',t,gmax*im_gr,'g-',t,gmax*ga,'b-',t,gmax*abs(gr),'k--')
plot(t,gmax*gx,'r-',t,gmax*gy,'g-',t,gmax*gz,'b-')
ylabel('g / Tm^-^1')
title(['b = ' num2str(b/1e9,2) '\cdot10^9 sm^-^2   q = ' num2str(qmax/2/pi/1e6,2) '\cdot10^6 m^-^1'])

subplot(2,2,2)
plot(t,dre_grdt,'r-',t,dim_grdt,'g-',t,dgadt,'b-',t,dgrdt,'k--')
xlabel('t / s'), ylabel('(dg/dt) / Tm^-^1s^-^1')

subplot(2,2,3)
plot(t,qx,'r-',t,qy,'g-',t,qz,'b-')
ylabel('q / m^-^1')
%title(['b = ' num2str(b/1e9,2) '\cdot10^9 sm^-^2   q = ' num2str(qmax/2/pi/1e6,2) '\cdot10^6 m^-^1'])

si = 16*1024;
Fx = fftshift(fft(qx,si));
Fy = fftshift(fft(qy,si));
Fz = fftshift(fft(qz,si));
Fx2 = abs(Fx).^2;
Fy2 = abs(Fy).^2;
Fz2 = abs(Fz).^2;
Fr2 = Fx2 + Fy2 + Fz2;

swh = 1/dt;
nu = swh*linspace(0,1,si); %frequency vector
nu = nu - nu(si/2+1); %set the proper zero frequency
subplot(2,2,4)
plot(nu,Fx2,'r-',nu,Fy2,'g-',nu,Fz2,'b-',nu,Fr2,'k--')
xlabel('\nu / s^-^1')
set(gca,'XLim',1000*[-1 1])

% Save waveforms in Bruker format
[out_path,out_name,out_ext] = fileparts(out_fn);
if ~isdir(out_path)
    mkdir(out_path)
end

title_str = ['Waveform for axisymmetric diffusion encoding. b_delta = ' num2str(b_delta,3)];

fid = fopen(out_fn,'w');

text.header = {
{['##TITLE= ' title_str]};
{'##JCAMP-DX= 5.00 Bruker JCAMP library'}
{'##DATA TYPE= Shape Data'}
{'##ORIGIN= Bruker Analytik GmbH'}
{'##OWNER= <nmrsu>'}
{['##DATE= ' datestr(now,'yyyy-mm-dd')]}
{['##TIME= ' datestr(now,'HH:MM:SS')]}
{'##MINX= 0'}
{'##MAXX= 1'}
{'##MINY= 0'}
{'##MAXY= 1'}
{'##$SHAPE_EXMODE= Gradient'}
{'##$SHAPE_TOTROT= 0'}
{'##$SHAPE_BWFAC= 0'}
{'##$SHAPE_INTEGFAC= 0'}
{'##$SHAPE_MODE= 0'}};

[nlines, ~] = size(text.header);

for nline = 1:nlines
    fprintf(fid,'%s\n',text.header{nline}{1});
end

fprintf(fid,'%s\n',['##BFACTOR=' num2str(bfactor)]);
fprintf(fid,'%s\n',['##DURATION=' num2str(tau)]);
fprintf(fid,'%s\n',['##DIRECTIONVEC= 0 0 1']);
fprintf(fid,'%s\n',['##NPOINTS=' num2str(numel(gx))]);
fprintf(fid,'%s\n',['##XYDATA= (T X Y Z)']);

formatspec = '%8.6f %8.6f %8.6f\r\n';

out_mat = [gx gy gz];
fprintf(fid, formatspec, out_mat');

fprintf(fid,'%s\n',['##END']);
fclose(fid);

max(sqrt(gx(:).^2+gy(:).^2+gz(:).^2))
