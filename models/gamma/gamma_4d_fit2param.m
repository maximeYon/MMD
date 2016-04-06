function res = gamma_4d_fit2param(mfs_fn, o_fn, opt)
% function fn = gamma_4d_fit2param(mfs_fn, o_path, opt)

if (nargin < 3), opt = []; end

res = -1;
    
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

% create parameter maps and save them

mfs.s0 = mfs.m(:,:,:,1);
mfs.iso = mfs.m(:,:,:,2);
mfs.mu2iso = mfs.m(:,:,:,3);
mfs.mu2aniso = mfs.m(:,:,:,4);

mfs.vlambda = 5/2*mfs.mu2aniso;
mfs.vlambda(isnan(mfs.vlambda)) = 0;
mfs.ufa = sqrt(3/2)*sqrt(1./(mfs.iso.^2./mfs.vlambda+1));
mfs.ufa(isnan(mfs.ufa)) = 0;
mfs.ciso = mfs.mu2iso./mfs.iso.^2;
mfs.ciso(isnan(mfs.ciso)) = 0;
mfs.cmu = mfs.ufa.^2;

mfs_fn = mdm_mfs_save(mfs, mfs.s, o_fn, opt);

res = 1;

