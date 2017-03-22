function dps = dtd_pake_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_pake_4d_fit2param(dps_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

    
opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them
dps.s0 = dps.m(:,:,:,1);
dps.iso = dps.m(:,:,:,2);
dps.delta = dps.m(:,:,:,3);

dps.par = dps.iso.*(1 + 2*dps.delta);
dps.perp = dps.iso.*(1 - dps.delta);
dps.aniso = dps.iso.*dps.delta;

dps.saniso = dps.aniso.^2;  % square anisotropic diffusivity
dps.saniso_n = dps.saniso./dps.iso.^2; % Normalized
dps.vlambda = 2*dps.saniso;
dps.ufa = sqrt(3/2)*sqrt(1./(dps.iso.^2./dps.vlambda+1));

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end



