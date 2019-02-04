function xps = mdm_bruker_dt_axderare2d_acqus2xps(data_path, xps_fn)
% function xps = mdm_bruker_dt_axderare2d_acqus2xps(data_path, xps_fn)
%
% Calculation of b-tensors for the Bruker pulse program DT_axderare2d
% as used in Topgaard, Phys. Chem. Chem. Phys., 2016.
% http://dx.doi.org/10.1039/c5cp07251d
%
% data_path: directory for gradient text files
% xps_fn (optional) filename of xps

if (nargin < 2), xps_fn = []; end

% Load Bruker acquisition parameters
load(fullfile(data_path,'NMRacqus'))

% xps calculation assumes that the pulse program is DT_axderare2d
if any(strcmp(NMRacqus.pulprog,{'DT_axderare2d','SM_moacrare2drs'})) ~= 1
    error('Wrong pulse program')
end

% Load gyromagnetic ratio gamma
gamma = mdm_bruker_gamma(NMRacqus);

% Load max gradient Gmax
Gmax = mdm_bruker_maxgradient(NMRacqus);

% Load number of images in indirect dimension td1
load(fullfile(data_path,'NMRacqu2s'))
td1 = NMRacqu2s.td;

% Index for powder averaging xps.a_ind
xps = mdm_xps_load(fullfile(data_path, 'g_xps.mat'));

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


% Read gradient ramps
% ramp.xa etc maps gradient modulations (a,b,c) to channels (x,y,z)
% Normalized from -1 to +1
gnams = {'xa','xb','xc','ya','yb','yc','za','zb','zc'};
for ngnam = 1:numel(gnams)
    gnam = gnams{ngnam};
    fn = fullfile(data_path,['g' gnam]);
    g = mdm_bruker_grad_read(fn);
    ramp.(gnam) = g;
end

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
G.x = repmat(G.a,[1, td1]).*repmat(G.xa',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.xb',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.xc',[Ndt, 1]);
G.y = repmat(G.a,[1, td1]).*repmat(G.ya',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.yb',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.yc',[Ndt, 1]);
G.z = repmat(G.a,[1, td1]).*repmat(G.za',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.zb',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.zc',[Ndt, 1]);

% Dephasing vector F.x in SI
F.x = cumsum(G.x*dt,1);
F.y = cumsum(G.y*dt,1);
F.z = cumsum(G.z*dt,1);
F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2); % Magnitude

% Diffusion weighting matrix b
% Factor 2 from the double DW blocks in DT_qVASrare2d
bt.xx = 2*gamma^2*sum(F.x.*F.x*dt,1)';
bt.xy = 2*gamma^2*sum(F.x.*F.y*dt,1)';
bt.xz = 2*gamma^2*sum(F.x.*F.z*dt,1)';
bt.yy = 2*gamma^2*sum(F.y.*F.y*dt,1)';
bt.yz = 2*gamma^2*sum(F.y.*F.z*dt,1)';
bt.zz = 2*gamma^2*sum(F.z.*F.z*dt,1)';

b = (bt.xx + bt.yy + bt.zz); % trace

% Timing parameters
tperacq = 15*NMRacqus.d22 + NMRacqus.vdT1 + 2e-6 + 20*NMRacqus.d32 + 4e-6*NMRacqus.p11 + ...
    NMRacqus.d34 + NMRacqus.vdT2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + ...
    NMRacqus.d35 + NMRacqus.aq + NMRacqus.d11;
indx = 1 + NMRacqus.nbl*(1:(NMRacqus.l2-1));
tsatperacq = tperacq - (.5e-6*NMRacqus.p11+3*NMRacqus.d32+2*NMRacqus.d22+NMRacqus.d11);
tperacq(indx) = tperacq(indx) + (NMRacqus.d22 + NMRacqus.d1);
t0acq = [0; cumsum(tperacq(1:(td1-1)))];
tT2W = .5e-6*NMRacqus.p11 + 7*NMRacqus.d32 + 2*NMRacqus.d22 + NMRacqus.vdT2 + NMRacqus.d34 + ...
    2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + NMRacqus.d35 + 0*NMRacqus.aq;
tT1R = zeros(td1,1);

Nslices = numel(unique(NMRacqus.fq1));
NRDenc = td1/Nslices;
[~,~,slice_ind] = unique(NMRacqus.fq1, 'rows');
for nslice = 1:Nslices
    indx = find(nslice == slice_ind);
    tT1R(indx(1)) = 1/eps;
    for n = 2:NRDenc
        tT1R(indx(n)) = t0acq(indx(n)) - (t0acq(indx(n-1))+tsatperacq(indx(n-1))) + 5*NMRacqus.d22 + NMRacqus.vdT1(indx(n)) + ...
            2e-6 + NMRacqus.d32 + .5e-6*NMRacqus.p11;
    end
end                 

% Save as xps
xps.b = b;
xps.n = td1;
xps.bt = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];
%xps.te = tT2W;
xps.te = NMRacqus.vdT2;

% xps.gmx = G.x'; % Maybe save full gradient modulations in the future
% xps.gmy = G.y';
% xps.gmz = G.z';

xps_old = xps;
xps = mdm_xps_calc_btpars(xps_old);

if any(strcmp(NMRacqus.pulprog,{'SM_moacrare2drs'}))
    ind = zeros(td1,1);
    ind(1:NMRacqus.nbl:(1+td1-NMRacqus.nbl)) = 1;
    xps = mdm_xps_subsample(xps, ind==1);
end
 

if (~isempty(xps_fn)), save(xps_fn,'xps'); end

