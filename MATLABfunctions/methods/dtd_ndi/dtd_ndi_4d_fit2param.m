function dps = dtd_ndi_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = kaden16_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

    
opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them
dps.s0     = dps.m(:,:,:,1);
dps.v_int  = dps.m(:,:,:,2);
dps.v_csf  = dps.m(:,:,:,3);

% suppress the v_int at high v_csf values -- it tends to get unreasonable
% there
if (opt.dtd_ndi.adjust_v_int)
    % with 0.4 and 0.1, we supress v_int when v_csf above 0.6
    dps.v_int = dps.v_int .* normcdf(1 - dps.v_csf, 0.4, 0.1);
end

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end

