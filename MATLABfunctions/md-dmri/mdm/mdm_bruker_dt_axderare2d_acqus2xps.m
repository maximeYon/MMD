function xps = mdm_bruker_dt_axderare2d_acqus2xps(data_path, xps_fn, rps)
% function xps = mdm_bruker_dt_axderare2d_acqus2xps(data_path, xps_fn)
%
% Calculation of b-tensors for the Bruker pulse program DT_axderare2d
% as used in Topgaard, Phys. Chem. Chem. Phys., 2016.
% http://dx.doi.org/10.1039/c5cp07251d
%
% data_path: directory for gradient text files
% xps_fn (optional) filename of xps

if (nargin < 2), xps_fn = []; end
if (nargin < 3), rps.maxomega = 500*2*pi; end

% Load Bruker acquisition parameters
load(fullfile(data_path,'NMRacqus'))

% xps calculation assumes that the pulse program is DT_axderare2d
if any(strcmp(NMRacqus.pulprog,{'DT_axderare2d','SM_moacrare2drs','DT_axderare2d_av4','DT_trtevcaxderare2d'})) ~= 1
    error('Wrong pulse program')
end

% Load gyromagnetic ratio gamma
gamma = mdm_bruker_gamma(NMRacqus);

% Load max gradient Gmax
Gmax = mdm_bruker_maxgradient(NMRacqus);

% Load number of images in indirect dimension td1
load(fullfile(data_path,'NMRacqu2s'))
td1 = NMRacqu2s.td;

% % Index for powder averaging xps.a_ind
% xps = mdm_xps_load(fullfile(data_path, 'g_xps.mat'));

% Read gradient time-modulations (a,b,c)
% Normalized from -1 to +1
gnams = {'a','b','c'};
for ngnam = 1:numel(gnams)
    gnam = gnams{ngnam};
    fn = fullfile(data_path,['g' gnam]);
    g = mdm_bruker_grad_read(fn);
    G.(gnam) = g;
end
Ndt = length(G.c); % Ndt: number of time steps
tau = NMRacqus.d3; % tau: duration of one full gradient modulation
dt = tau/Ndt; % dt: time step
G_pre = G;

Ndt_d4 = round(NMRacqus.d4/dt);
Ndt_d32 = round(NMRacqus.d32/dt);
Ndt_p12 = round(NMRacqus.p12/1e6/dt); Ndt_p12 = 2*ceil(Ndt_p12/2);

ind_slice = logical([zeros(Ndt,1); ones(Ndt_d4+4*Ndt_d32+Ndt_p12,1); zeros(Ndt,1)]);
Ndt0 = Ndt;
Ndt = numel(ind_slice);
ind_pre180 = logical([ones(Ndt0+Ndt_d4+2*Ndt_d32+Ndt_p12/2,1); zeros(Ndt0+2*Ndt_d32+Ndt_p12/2,1)]); 
ind_pre180 = repmat(ind_pre180,[1, td1]);

Gslice.x = [zeros(Ndt_d4,1); Gmax.*NMRacqus.cnst34/100*linspace(0,1,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst34/100*linspace(1,0,Ndt_d32)' + Gmax.*NMRacqus.cnst31/100*linspace(0,1,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst31/100*ones(Ndt_p12,1); ...
     Gmax.*NMRacqus.cnst34/100*linspace(0,1,Ndt_d32)' + Gmax.*NMRacqus.cnst31/100*linspace(1,0,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst34/100*linspace(1,0,Ndt_d32)'];        
Gslice.y = [zeros(Ndt_d4,1); Gmax.*NMRacqus.cnst35/100*linspace(0,1,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst35/100*linspace(1,0,Ndt_d32)' + Gmax.*NMRacqus.cnst32/100*linspace(0,1,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst32/100*ones(Ndt_p12,1); ...
     Gmax.*NMRacqus.cnst35/100*linspace(0,1,Ndt_d32)' + Gmax.*NMRacqus.cnst32/100*linspace(1,0,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst35/100*linspace(1,0,Ndt_d32)'];        
Gslice.z = [zeros(Ndt_d4,1); Gmax.*NMRacqus.cnst36/100*linspace(0,1,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst36/100*linspace(1,0,Ndt_d32)' + Gmax.*NMRacqus.cnst33/100*linspace(0,1,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst33/100*ones(Ndt_p12,1); ...
     Gmax.*NMRacqus.cnst36/100*linspace(0,1,Ndt_d32)' + Gmax.*NMRacqus.cnst33/100*linspace(1,0,Ndt_d32)'; ...
     Gmax.*NMRacqus.cnst36/100*linspace(1,0,Ndt_d32)'];        

G.a = [G.a; zeros(Ndt_d4+4*Ndt_d32+Ndt_p12,1); G.a];
G.b = [G.b; zeros(Ndt_d4+4*Ndt_d32+Ndt_p12,1); G.b];
G.c = [G.c; zeros(Ndt_d4+4*Ndt_d32+Ndt_p12,1); G.c];

% Read gradient ramps
% ramp.xa etc maps gradient modulations (a,b,c) to channels (x,y,z)
% Normalized from -1 to +1
gnams = {'xa','xb','xc','ya','yb','yc','za','zb','zc'};
for ngnam = 1:numel(gnams)
    gnam = gnams{ngnam};
    fn = fullfile(data_path,['g' gnam]);
    g = mdm_bruker_grad_read(fn);
    ramp.(gnam) = g(1:td1);
end

% Convert ramps r.ax to uvec axis of q-vector cone
gpas.avecnorm = sqrt(ramp.xa.^2+ramp.ya.^2+ramp.za.^2);
gpas.bvecnorm = sqrt(ramp.xb.^2+ramp.yb.^2+ramp.zb.^2);
gpas.cvecnorm = sqrt(ramp.xc.^2+ramp.yc.^2+ramp.zc.^2);
gpas.avec = [ramp.xa ramp.ya ramp.za]./repmat(gpas.avecnorm,[1 3]);
gpas.bvec = [ramp.xb ramp.yb ramp.zb]./repmat(gpas.bvecnorm,[1 3]);
gpas.cvec = [ramp.xc ramp.yc ramp.zc]./repmat(gpas.cvecnorm,[1 3]);
gpas.abcrossvec = msf_crossprod_nx3vectors(gpas.avec,gpas.bvec);
uvec = zeros(td1,3);
ind = find(gpas.avecnorm>0);
uvec(ind,:) = gpas.abcrossvec(ind,:);
ind = find(gpas.cvecnorm>0);
uvec(ind,:) = gpas.cvec(ind,:);

% Convert normlized ramps r.ax to G.xa in SI 
G.xa = ramp.xa*Gmax.*NMRacqus.cnst1/100; % NMRacqus.cnst1: scaling factor for x-gradient. Max 100%. 
G.xb = ramp.xb*Gmax.*NMRacqus.cnst1/100;
G.xc = ramp.xc*Gmax.*NMRacqus.cnst1/100;
G.ya = ramp.ya*Gmax.*NMRacqus.cnst2/100; % NMRacqus.cnst2: scaling factor for y-gradient. Max 100%. 
G.yb = ramp.yb*Gmax.*NMRacqus.cnst2/100;
G.yc = ramp.yc*Gmax.*NMRacqus.cnst2/100;
G.za = ramp.za*Gmax.*NMRacqus.cnst3/100; % NMRacqus.cnst3: scaling factor for z-gradient. Max 100%. 
G.zb = ramp.zb*Gmax.*NMRacqus.cnst3/100;
G.zc = ramp.zc*Gmax.*NMRacqus.cnst3/100;

% Transformation from gradients (a,b,c) to (x,y,z)
% G.x: Ndt x td1 array in SI
%G.a = G_pre.a; G.b = G_pre.b; G.c = G_pre.c; Ndt = length(G.c);

G.x = repmat(G.a,[1, td1]).*repmat(G.xa',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.xb',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.xc',[Ndt, 1]);
G.y = repmat(G.a,[1, td1]).*repmat(G.ya',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.yb',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.yc',[Ndt, 1]);
G.z = repmat(G.a,[1, td1]).*repmat(G.za',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.zb',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.zc',[Ndt, 1]);

G.x(ind_slice,:) = repmat(Gslice.x,[1, td1]);
G.y(ind_slice,:) = repmat(Gslice.y,[1, td1]);
G.z(ind_slice,:) = repmat(Gslice.z,[1, td1]);

G.x(ind_pre180) = -G.x(ind_pre180);
G.y(ind_pre180) = -G.y(ind_pre180);
G.z(ind_pre180) = -G.z(ind_pre180);

% figure(1), clf, plot(1:Ndt,G.x,'r-',1:Ndt,G.y,'g-',1:Ndt,G.z,'b-'), pause(1)

% Dephasing vector F.x in SI
F.x = cumsum(G.x*dt,1);
F.y = cumsum(G.y*dt,1);
F.z = cumsum(G.z*dt,1);
F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2); % Magnitude

% Diffusion weighting matrix b
% Factor 2 from the double DW blocks in DT_axderare2d
% N_Gblocks = 2;
N_Gblocks = 1;
bt.xx = N_Gblocks*gamma^2*sum(F.x.*F.x*dt,1)';
bt.xy = N_Gblocks*gamma^2*sum(F.x.*F.y*dt,1)';
bt.xz = N_Gblocks*gamma^2*sum(F.x.*F.z*dt,1)';
bt.yy = N_Gblocks*gamma^2*sum(F.y.*F.y*dt,1)';
bt.yz = N_Gblocks*gamma^2*sum(F.y.*F.z*dt,1)';
bt.zz = N_Gblocks*gamma^2*sum(F.z.*F.z*dt,1)';

b = (bt.xx + bt.yy + bt.zz); % trace

% % Timing parameters
% tperacq = 15*NMRacqus.d22 + NMRacqus.d1 + 2e-6 + 20*NMRacqus.d32 + 1e-6*NMRacqus.p11 + ...
%     NMRacqus.d34 + NMRacqus.vdT2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + ...
%     NMRacqus.d35 + NMRacqus.aq + NMRacqus.d11;
% indx = 1 + NMRacqus.nbl*(1:(NMRacqus.l2-1));
% tsatperacq = tperacq - (.5e-6*NMRacqus.p11+3*NMRacqus.d32+2*NMRacqus.d22+NMRacqus.d11);
% tperacq(indx) = tperacq(indx) + (NMRacqus.d22 + NMRacqus.d1);
% t0acq = [0; cumsum(tperacq(1:(td1-1)))];
% tT2W = .5e-6*NMRacqus.p11 + 7*NMRacqus.d32 + 2*NMRacqus.d22 + NMRacqus.vdT2 + NMRacqus.d34 + ...
%     2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + NMRacqus.d35 + 0*NMRacqus.aq;
% tT1R = zeros(td1,1);

NMRacqus.dwcomp = NMRacqus.dw - 1.25e-6;
tperloop = 6*NMRacqus.dw + 2*NMRacqus.dwcomp + 1e-6*NMRacqus.p12 + 4*NMRacqus.d32 + NMRacqus.d40;
tT2W = .5e-6*NMRacqus.p11 + 10*NMRacqus.d32 + NMRacqus.d34 + 2*NMRacqus.d22 + 2*NMRacqus.d6  + ...
    2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + 4*NMRacqus.dw + NMRacqus.d41 + ...
    .5*NMRacqus.l1*tperloop + tperloop - 2*NMRacqus.d32 - .5*NMRacqus.d40;
tT1R = .5*1e-6*NMRacqus.p11 + 3*10*NMRacqus.d32 + 2*NMRacqus.d22 + NMRacqus.d11 + NMRacqus.d1 + 2e-3;


% Save as xps
xps.b = b;
xps.n = td1;
xps.bt = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];
xps.u = uvec;
xps.theta = acos(uvec(:,3));
xps.phi = atan2(uvec(:,2),uvec(:,1));
xps.te = tT2W*ones(xps.n,1);
xps.tr = tT1R*ones(xps.n,1);

if any(strcmp(NMRacqus.pulprog,{'DT_trtevcaxderare2d'})) == 1
    xps.te = .5e-6*NMRacqus.p11 + 10*NMRacqus.d32 + NMRacqus.d34 + 2*NMRacqus.d22 + 2*.5*NMRacqus.vdte(1:td1)  + ...
        NMRacqus.vc(1:td1)*2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + 4*NMRacqus.dw + NMRacqus.d41 + ...
        .5*NMRacqus.l1*tperloop + tperloop - 2*NMRacqus.d32 - .5*NMRacqus.d40;
    xps.tr = .5*1e-6*NMRacqus.p11 + 3*NMRacqus.d32 + 2*NMRacqus.d22 + NMRacqus.d11 + ...
        2e-3 + NMRacqus.vdtr(1:td1) + NMRacqus.d22 + 1e-6*NMRacqus.de + 0.25e-6 + NMRacqus.d32 + .5*1e-6*NMRacqus.p11;
end

xps_temp = xps;

xps_temp.gwf_x = G.x'; % Maybe save full gradient modulations in the future
xps_temp.gwf_y = G.y';
xps_temp.gwf_z = G.z';
xps_temp.gwf_t = repmat(dt*(1:1:Ndt),[td1 1]);
%figure(1), clf, plot(xps_temp.gwf_t,xps_temp.gwf_x,'r-')

xps_temp = mdm_xps_calc_btpars(xps_temp);
xps_temp = mdm_xps_add_btomega(xps_temp,rps);

xps.b_delta = xps_temp.b_delta;
xps.b_eta = xps_temp.b_eta;
xps.btomega = xps_temp.btomega;
xps.domega = xps_temp.domega;
xps.rmsomega = xps_temp.rmsomega;
xps.momega = xps_temp.momega;

figure(11), clf, plot((0:(size(xps.btomega,2)-1))',xps.btomega','-'), pause(.1)

if (~isempty(xps_fn)), save(xps_fn,'xps'); end

