function dps = dtd_covariance_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_covariance_4d_fit2param(mfs_fn, dps_fn, opt)
%
% Compute dps (display parameters structure) from the mfs (model fit
% structure) using functions tm_dt_to_dps and tm_ct_to_dps

if (nargin < 2), mfs_fn = []; end
if (nargin < 3), opt    = []; end

opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
sz  = msf_size(mfs.m(:,:,:,1), 3);

% reshape help functions
g = @(a,n) reshape(a, prod(sz(1:3)), n);
f = @(a,n) reshape(a, sz(1), sz(2), sz(3), n);

% init dps
dps.nii_h = mfs.nii_h;
dps.mask  = mfs.mask;

% get data from input
dps.s0  = mfs.m(:,:,:,1);

% diffusion tensor
dt_1x6 = g(mfs.m(:,:,:,2:7), 6) * 1e9;
dps   = tm_dt_to_dps(dt_1x6, dps, f);

% covariance tensor
ct_1x21  = g(mfs.m(:,:,:,8:28), 21) * 1e18;
dps = tm_ct_to_dps(ct_1x21, dps, f);

if (~isempty(dps_fn))
    mdm_dps_save(dps, mfs.s, dps_fn, opt);
end


