function res = dti_euler_4d_fit2param(mfs_fn, o_fn, opt)
% function fn = dti_euler_4d_fit2param(mfs_fn, o_path, opt)

if (nargin < 3), opt = []; end

res = -1;
    
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

% create parameter maps and save them

mfs.s0 = mfs.m(:,:,:,1);
mfs.lambdax = mfs.m(:,:,:,2);
mfs.lambday = mfs.m(:,:,:,3);
mfs.lambdaz = mfs.m(:,:,:,4);
mfs.euler_alpha = angle(exp(i*mfs.m(:,:,:,5)));
mfs.euler_beta = angle(exp(i*mfs.m(:,:,:,6)));
mfs.euler_gamma = angle(exp(i*mfs.m(:,:,:,7)));

sz = [1 1 1];
sz_temp = size(mfs.s0);
sz(1:numel(sz_temp)) = sz_temp;
mfs.t1x6 = zeros([sz(1) sz(2) sz(3) 6]);
mfs.lambdazzvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambdaxxvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambdayyvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambda11vec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambda22vec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambda33vec = zeros([sz(1) sz(2) sz(3) 3]);
param = {'trace','iso','lambda33','lambda22','lambda11','lambdazz','lambdaxx','lambdayy','vlambda',...
    'delta','eta','s','p','l','fa','cs','cl','cp','cm'};
for nparam = 1:numel(param)
    eval(['mfs.' param{nparam} ' = zeros([sz(1) sz(2) sz(3)]);']);
end
for nk = 1:sz(3)
    for nj = 1:sz(2)
        for ni = 1:sz(1)
            alpha = mfs.euler_alpha(ni,nj,nk);
            beta = mfs.euler_beta(ni,nj,nk);
            gamma = mfs.euler_gamma(ni,nj,nk);
            [rotmat,rotmatinv] = tm_euler_angles2rotmat(alpha,beta,gamma);
            lambdax = mfs.lambdax(ni,nj,nk);
            lambday = mfs.lambday(ni,nj,nk);
            lambdaz = mfs.lambdaz(ni,nj,nk);
            dt_lambda = [
                lambdax 0 0
                0 lambday 0
                0 0 lambdaz];
            dt3x3 = rotmat*dt_lambda*rotmatinv;
            dt = tm_t2tpars(dt3x3);
            mfs.t1x6(ni,nj,nk,:) = dt.t1x6;
            mfs.lambdazzvec(ni,nj,nk,:) = dt.lambdazzvec;
            mfs.lambdaxxvec(ni,nj,nk,:) = dt.lambdaxxvec;
            mfs.lambdayyvec(ni,nj,nk,:) = dt.lambdayyvec;
            mfs.lambda11vec(ni,nj,nk,:) = dt.lambda11vec;
            mfs.lambda22vec(ni,nj,nk,:) = dt.lambda22vec;
            mfs.lambda33vec(ni,nj,nk,:) = dt.lambda33vec;
            for nparam = 1:numel(param)
                eval(['mfs.' param{nparam} '(ni,nj,nk) = dt.' param{nparam} ';']);
            end

        end
    end
end

mfs_fn = mdm_mfs_save(mfs, mfs.s, o_fn, opt);

res = 1;

