function res = mdm_bruker_axsymrare2d_acqus2xps(paths)
% function res = mdm_bruker_axsymrare2d_acqus2xps(paths)
%
% Calculation of b-tensors for the Bruker pulse program DT_qVASrare2d
% as used in Topgaard, Phys. Chem. Chem. Phys., 2016.
% http://dx.doi.org/10.1039/c5cp07251d
%
% paths.data: directory for gradient text files
% paths.out: directory for experimental parameter structure xps.mat
%

% Load Bruker acquisition parameters
load(fullfile(paths.data,'NMRacqus'))

% xps calculation assumes that the pulse program is DT_qVASrare2d
if any(strcmp(NMRacqus.pulprog,{'DT_qVASrare2d'})) ~= 1
    error('Wrong pulse program')
end

% Load gyromagnetic ratio gamma
gamma = mdm_bruker_gamma(NMRacqus);

% Load max gradient Gmax
Gmax = mdm_bruker_maxgradient(NMRacqus);

% Load number of images in indirect dimension td1
load(fullfile(paths.data,'NMRacqu2s'))
td1 = NMRacqu2s.td;

% Index for powder averaging xps.a_ind
load(fullfile(paths.data,'/DiffRamp'))
[b_ind,bd_ind,bu_ind] = ndgrid(1:DiffRamp.NbTtrace,1:DiffRamp.NbTDelta,1:DiffRamp.NbTdir);
xps.b_ind = reshape(b_ind,td1,1); % Same b_ind means same b-tensor size
xps.bd_ind = reshape(bd_ind,td1,1); % Same bd_ind means same b-tensor shape
xps.bu_ind = reshape(bu_ind,td1,1); % Same bu_ind means same b-tensor orientation
[a_ind,~] = ndgrid(1:(DiffRamp.NbTtrace*DiffRamp.NbTDelta),1:DiffRamp.NbTdir);
xps.a_ind = reshape(a_ind,td1,1); % Same a_ind means same b-tensor size and shape

% Read qVAS gradient time-modulations (a,b,c)
% Normalized from -1 to +1
gnams = {'a','b','c'};
for ngnam = 1:numel(gnams)
    gnam = gnams{ngnam};
    fn = fullfile(paths.data,['qVAS' gnam]);
    g = mdm_bruker_readgrad(fn);
    eval(['G.' gnam ' = g;'])
end
Ndt = length(G.c); % Ndt: number of time steps
tau = NMRacqus.d3; % tau: duration of one full gradient modulation
dt = tau/Ndt; % dt: time step
%figure(1), clf, plot((1:Ndt)',[G.c G.b G.a],'-'), return

% Read gradient ramps
% ramp.ax etc maps gradient modulations (a,b,c) to channels (x,y,z)
% Normalized from -1 to +1
gnams = {'ax','ay','az','bx','by','bz','cx','cy','cz'};
for ngnam = 1:numel(gnams)
    gnam = gnams{ngnam};
    fn = fullfile(paths.data,['r' gnam]);
    g = mdm_bruker_readgrad(fn);
    eval(['ramp.' gnam ' = g;'])
end

% Convert normlized ramps r.ax to G.ax in SI 
G.ax = ramp.ax*Gmax.*NMRacqus.cnst1/100; % NMRacqus.cnst1: scaling factor for x-gradient. Max 100%. 
G.bx = ramp.bx*Gmax.*NMRacqus.cnst1/100;
G.cx = ramp.cx*Gmax.*NMRacqus.cnst1/100;
G.ay = ramp.ay*Gmax.*NMRacqus.cnst2/100; % NMRacqus.cnst2: scaling factor for y-gradient. Max 100%. 
G.by = ramp.by*Gmax.*NMRacqus.cnst2/100;
G.cy = ramp.cy*Gmax.*NMRacqus.cnst2/100;
G.az = ramp.az*Gmax.*NMRacqus.cnst3/100; % NMRacqus.cnst3: scaling factor for z-gradient. Max 100%. 
G.bz = ramp.bz*Gmax.*NMRacqus.cnst3/100;
G.cz = ramp.cz*Gmax.*NMRacqus.cnst3/100;

% Transformation from gradients (a,b,c) to (x,y,z)
% G.x: Ndt x td1 array in SI
G.x = repmat(G.a,[1, td1]).*repmat(G.ax',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.bx',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.cx',[Ndt, 1]);
G.y = repmat(G.a,[1, td1]).*repmat(G.ay',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.by',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.cy',[Ndt, 1]);
G.z = repmat(G.a,[1, td1]).*repmat(G.az',[Ndt, 1]) + ...
    repmat(G.b,[1, td1]).*repmat(G.bz',[Ndt, 1]) + ...
    repmat(G.c,[1, td1]).*repmat(G.cz',[Ndt, 1]);

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

% Save as xps
xps.b = b;
xps.n = td1;
xps.bt = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];
% xps.gmx = G.x'; % Maybe save full gradient modulations in the future
% xps.gmy = G.y';
% xps.gmz = G.z';

xps_old = xps;
xps = mdm_xps_bt2btpars(xps_old);

save(paths.xps_fn,'xps')

res = 1;
