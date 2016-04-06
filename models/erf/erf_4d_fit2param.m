function res = erf_4d_fit2param(mfs_fn, o_fn, opt)
% function fn = erf_4d_fit2param(mfs_fn, o_path, opt)

if (nargin < 3), opt = []; end

res = -1;
    
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

% create parameter maps and save them

mfs.s0 = mfs.m(:,:,:,1);
mfs.iso = mfs.m(:,:,:,2);
mfs.delta = mfs.m(:,:,:,3);

mfs.par = mfs.iso.*(1 + 2*mfs.delta);
mfs.perp = mfs.iso.*(1 - mfs.delta);
mfs.logratio = log10(mfs.par./mfs.perp);
mfs.logratio(isnan(mfs.logratio)) = 0;
mfs.vlambda = 2*(mfs.iso.*mfs.delta).^2;
mfs.vlambda(isnan(mfs.vlambda)) = 0;
mfs.mu2aniso = 2/5*mfs.vlambda;
mfs.ufa = sqrt(3/2)*sqrt(1./(mfs.iso.^2./mfs.vlambda+1));
mfs.ufa(isnan(mfs.ufa)) = 0;
mfs.cmu = mfs.ufa.^2;

mfs_fn = mdm_mfs_save(mfs, mfs.s, o_fn, opt);

res = 1;

