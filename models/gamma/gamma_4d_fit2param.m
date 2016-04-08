function dps = gamma_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = gamma_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them
dps.s0 = dps.m(:,:,:,1);
dps.iso = dps.m(:,:,:,2);
dps.mu2iso = dps.m(:,:,:,3);
dps.mu2aniso = dps.m(:,:,:,4);

dps.vlambda = 5/2*dps.mu2aniso;
dps.vlambda(isnan(dps.vlambda)) = 0;
dps.ufa = sqrt(3/2)*sqrt(1./(dps.iso.^2./dps.vlambda+1));
dps.ufa(isnan(dps.ufa)) = 0;
dps.ciso = dps.mu2iso./dps.iso.^2;
dps.ciso(isnan(dps.ciso)) = 0;
dps.cmu = dps.ufa.^2;

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end
