function xps = mdm_xps_add_btomega(xps_in,rps)
% Add tensor-valued encoding spectrum b(omega) to xps
%
% Input
% xps_in: eXperimental Parameter Structure with fields gwf_x, gwf_y, gwf_z, gwf_t
% with effective (including 180 pulse reversals) gradient waveforms (x,y,z) and time vector,
% all as arrays with dimensions [xps.n, number of time steps]
% Also requires input field b with b-value from conventional time-domain
% calculation.
% rps: Recon Parameter Structure with fields maxomega and nomega for
% generating output omega vector by maxomega*linspace(0,1,nomega)
%
% Output
% xps with added fields 
% btomega: tensor-valued encoding spectrum in pseudo-Voigt format with dimensions [xps.n, 6*nomega]
% domega: spectral resolution in omega dimension [xps.n, 1]
% momega: centroid encoding frequency (approximately equal to
% main encoding peaks for q-MAS and q-DOR) [xps.n, 1]

xps = xps_in;

g.x = xps.gwf_x;
g.y = xps.gwf_y;
g.z = xps.gwf_z;
t = xps.gwf_t(1,:);

Nt = size(g.x,2); % Number of time steps in waveform
Nenc = size(g.x,1); % Number of encodings
dt = t(1,2) - t(1,1); %Time steps. Assumes same t vector for all encodings

% Load gyromagnetic ratio gamma
NMRacqus.nuc1 = '1H';
gamma = mdm_bruker_gamma(NMRacqus);

% Dephasing vector q.x in SI
q.x = gamma*cumsum(g.x*dt,2);
q.y = gamma*cumsum(g.y*dt,2);
q.z = gamma*cumsum(g.z*dt,2);

% Diffusion encoding spectrum des
% First temporary wide-range and high-resolution spectrum
domega_interp = rps.maxomega/rps.nomega;
maxt_interp = 2*pi/domega_interp;
Nt_interp = maxt_interp/dt;
Nomega = max([Nt pow2(nextpow2(Nt_interp)+1)]);

omega = 2*pi/dt*linspace(0,1,Nomega); omega = omega - omega(Nomega/2+1);
domega = omega(2)-omega(1);
qomega.x = fftshift(fft(q.x,Nomega,2),2)*dt; 
qomega.y = fftshift(fft(q.y,Nomega,2),2)*dt; 
qomega.z = fftshift(fft(q.z,Nomega,2),2)*dt; 
des.sxx = qomega.x.*(conj(qomega.x));
des.syy = qomega.y.*(conj(qomega.y));
des.szz = qomega.z.*(conj(qomega.z));
des.sxy = qomega.x.*(conj(qomega.y));
des.sxz = qomega.x.*(conj(qomega.z));
des.syz = qomega.y.*(conj(qomega.z));
des.strace = des.sxx + des.syy + des.szz;

des.omega = omega;
des.domega = domega;

rmsomega = sqrt(sum(repmat(des.omega,[size(des.strace,1) 1]).^2.*des.strace,2)./sum(des.strace,2));

omegalim = rps.maxomega;
if omegalim < 1.5*max(rmsomega(:))
    warning(['rps.maxomega/2pi = ' num2str(omegalim/2/pi) ' Hz should be much larger than rmsomega/2pi = ' num2str(max(rmsomega(:))/2/pi) ' Hz'])
end

% Truncate at max frequency and remove redundant values at negative
% frequencies
ind = (omega < omegalim) & (omega > -eps); 
fields = {'sxx','syy','szz','sxy','sxz','syz','strace','omega'};
for nfield = 1:numel(fields)
    field = fields{nfield};
    des.(field) = des.(field)(:,ind);
end

% Downsample to the omega vector requested by input rps
des_interp.omega = rps.maxomega*linspace(0,1,rps.nomega);
fields = {'sxx','syy','szz','sxy','sxz','syz'};
for nfield = 1:numel(fields)
    field = fields{nfield};
    des_interp.(field) = interp1(des.omega',des.(field)',des_interp.omega','makima')';
end
des_interp.strace = des_interp.sxx+des_interp.syy+des_interp.szz;

%Uncomment to check quality of downsampling
% size(des.strace), size(des_interp.strace)
% figure(1), clf, plot(des.omega'/2/pi,des.strace','-',des_interp.omega'/2/pi,des_interp.strace','o'), return

des = des_interp;

% Multiply zero-value with 0.5. Belongs half to positive and half to
% negative as in FFTs.
for nfield = 1:numel(fields)
    field = fields{nfield};
    des.(field)(:,1) = .5*des.(field)(:,1);
end

% Assemble to pseudo-Voigt format. Renormalize to eliminitate round-off
% errors in b-value
btomega = [des.sxx des.syy des.szz sqrt(2)*[des.sxy des.sxz des.syz]]...
    .*repmat(xps.b./sum(des.sxx+des.syy+des.szz,2),[1 6*rps.nomega]);

% Save relevant fields to output xps
xps.btomega = msf_notfinite2zero(real(btomega));
xps.domega = (des.omega(1,2)-des.omega(1,1))*ones(Nenc,1);
xps.rmsomega = msf_notfinite2zero(rmsomega);

Nomega = size(xps.btomega,2)/6;
omega = repmat(xps.domega,[1 Nomega]).*repmat(0:(Nomega-1),[xps.n 1]);
bomega = xps.btomega(:,1:Nomega) + xps.btomega(:,(1:Nomega)+Nomega) + xps.btomega(:,(1:Nomega)+2*Nomega);
xps.momega = double(sum(omega.*bomega,2)./sum(bomega,2));

