function res = saupe_4d_fit2param(paths, o_fn, opt)
% function fn = saupe_4d_fit2param(mfs_fn, o_path, opt)

if (nargin < 3), opt = []; end

res = -1;

opt = mdm_opt(opt);

% create parameter maps and save them

mfs_dti = mdm_mfs_load(paths.mfs.dti_derived_fn);
mfs_gamma = mdm_mfs_load(paths.mfs.gamma_derived_fn);
mfs_erf = mdm_mfs_load(paths.mfs.erf_derived_fn);

mfs = mdm_mfs_load(paths.mfs.dti_primary_fn);
sz = size(mfs.m);

kronecker = permute(repmat([1 1 1 0 0 0]',[1 sz(3) sz(2) sz(1)]),[4 3 2 1]);
mfs.t1x6 = (mfs_dti.t1x6./repmat(mfs_dti.iso,[1 1 1 6]) - kronecker)./repmat(mfs_erf.delta,[1 1 1 6])/2;
mfs.t1x6prim = (2*mfs.t1x6 + kronecker)/3;
mfs.t1x6(isnan(mfs.t1x6)) = 0;
mfs.t1x6prim(isnan(mfs.t1x6prim)) = 0;

mfs.slambdaxx = (mfs_dti.lambdaxx./mfs_dti.iso - 1)./mfs_erf.delta/2;
mfs.slambdaxx(isnan(mfs.slambdaxx)) = 0;
mfs.slambdayy = (mfs_dti.lambdayy./mfs_dti.iso - 1)./mfs_erf.delta/2;
mfs.slambdayy(isnan(mfs.slambdayy)) = 0;
mfs.slambdazz = (mfs_dti.lambdazz./mfs_dti.iso - 1)./mfs_erf.delta/2;
mfs.slambdazz(isnan(mfs.slambdazz)) = 0;

mfs.slambdaxxprim = (2*mfs.slambdaxx + 1)/3;
mfs.slambdayyprim = (2*mfs.slambdayy + 1)/3;
mfs.slambdazzprim = (2*mfs.slambdazz + 1)/3;

mfs.slambdas = zeros([sz(1) sz(2) sz(3) 3]);
mfs.slambdas(:,:,:,1) = mfs.slambdaxxprim;
mfs.slambdas(:,:,:,2) = mfs.slambdayyprim;
mfs.slambdas(:,:,:,3) = mfs.slambdazzprim;
mfs.slambda11prim = min(mfs.slambdas,[],4);
mfs.slambda33prim = max(mfs.slambdas,[],4);
mfs = rmfield(mfs,'slambdas');

mfs.cc = mfs_dti.cm./mfs_gamma.cmu;
mfs.cc(isnan(mfs.cc)) = 0;
mfs.cc(mfs.cc>1) = 1;
mfs.op = sqrt(mfs.cc);

mfs_fn = mdm_mfs_save(mfs, mfs.s, o_fn, opt);

res = 1;

