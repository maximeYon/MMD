function dps = dtd_ndi_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = kaden16_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

    
opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them
dps.s0     = dps.m(:,:,:,1);
dps.v_int  = dps.m(:,:,:,2);
dps.lambda = dps.m(:,:,:,3);

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end

