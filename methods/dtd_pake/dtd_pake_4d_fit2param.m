function dps = dtd_pake_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_pake_4d_fit2param(dps_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

    
opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them
dps.s0      = dps.m(:,:,:,1);
dps.iso     = dps.m(:,:,:,2) * 1e9;
dps.delta   = dps.m(:,:,:,3);

dps.par             = dps.iso.*(1 + 2*dps.delta);
dps.perp            = dps.iso.*(1 - dps.delta);
dps.logratio        = log10(dps.par./dps.perp);


dps.vlambda         = 2*(dps.iso.*dps.delta).^2;
dps.ufa             = sqrt(3/2)*sqrt(1./(dps.iso.^2./dps.vlambda+1));
dps.cmu             = dps.ufa.^2;

dps.vlambda(isnan(dps.vlambda)) = 0;
dps.logratio(isnan(dps.logratio)) = 0;
dps.ufa(isnan(dps.ufa)) = 0;


dps.mu2aniso        = 2/5*dps.vlambda;


if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end



