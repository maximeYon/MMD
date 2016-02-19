function fn = dti_euler_4d_fit2param(mfs_fn, o_path, opt)
% function fn = dti_euler_4d_fit2param(mfs_fn, o_path, opt)

if (nargin < 3), opt = []; end
    
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
mfs.dt.t1x6 = zeros([sz(1) sz(2) sz(3) 6]);
mfs.dt.lambdazzvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.dt.lambdaxxvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.dt.lambdayyvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.dt.lambdaminvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.dt.lambdamidvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.dt.lambdamaxvec = zeros([sz(1) sz(2) sz(3) 3]);

param = {'trace','iso','lambdamax','lambdamid','lambdamin','lambdazz','lambdaxx','lambdayy','vlambda',...
    'delta','eta','s','p','l','fa','cl','cp'};
for nparam = 1:numel(param)
    eval(['mfs.dt.' param{nparam} ' = zeros([sz(1) sz(2) sz(3)]);']);
end

for nk = 1:sz(3)
    for nj = 1:sz(2)
        for ni = 1:sz(1)
            alpha = mfs.euler_alpha(ni,nj,nk);
            beta = mfs.euler_beta(ni,nj,nk);
            gamma = mfs.euler_gamma(ni,nj,nk);
            [rotmat,rotmatinv] = euler_angles2rotmat(alpha,beta,gamma);
            lambdax = mfs.lambdax(ni,nj,nk);
            lambday = mfs.lambday(ni,nj,nk);
            lambdaz = mfs.lambdaz(ni,nj,nk);
            dt_lambda = [
                lambdax 0 0
                0 lambday 0
                0 0 lambdaz];
            dt3x3 = rotmat*dt_lambda*rotmatinv;
            dt = t2tpars(dt3x3);
            mfs.dt.t1x6(ni,nj,nk,:) = dt.t1x6;
            mfs.dt.lambdazzvec(ni,nj,nk,:) = dt.lambdazzvec;
            mfs.dt.lambdaxxvec(ni,nj,nk,:) = dt.lambdaxxvec;
            mfs.dt.lambdayyvec(ni,nj,nk,:) = dt.lambdayyvec;
            mfs.dt.lambdaminvec(ni,nj,nk,:) = dt.lambdaminvec;
            mfs.dt.lambdamidvec(ni,nj,nk,:) = dt.lambdamidvec;
            mfs.dt.lambdamaxvec(ni,nj,nk,:) = dt.lambdamaxvec;
            for nparam = 1:numel(param)
                eval(['mfs.dt.' param{nparam} '(ni,nj,nk) = dt.' param{nparam} ';']);
            end
        end
    end
end

n_map = 8;
fn = cell(1,n_map);
for c = 1:n_map
    
    min_max = [-inf inf];
    switch (c)
        
        case 1
            param = 's0';
            min_max = [0 inf];
            x = mfs.m(:,:,:,c);

        case 2
            param = 'lambdax';
            min_max = [0 10e-9];
            x = mfs.m(:,:,:,c);
            
        case 3
            param = 'lambday';
            min_max = [0 10e-9];
            x = mfs.m(:,:,:,c);
            
        case 4
            param = 'lambdaz';
            min_max = [0 10e-9];
            x = mfs.m(:,:,:,c);

        case 5
            param = 'euler_alpha';
            x = angle(exp(i*mfs.m(:,:,:,c)));

        case 6
            param = 'euler_beta';
            x = angle(exp(i*mfs.m(:,:,:,c)));
            
        case 7
            param = 'euler_gamma';
            x = angle(exp(i*mfs.m(:,:,:,c)));
            
        case 8
            param = 'md';
            min_max = [0 10e-9];
            x = sum(mfs.m(:,:,:,2:4),4);
    end
    
    fn{c} = fullfile(o_path, ['dti_euler_' param opt.nii_ext]);
    
%     % make sure the min_max field is there
%     opt.vasco16.(param).present = 1;
%     opt.vasco16.(param) = msf_ensure_field(opt.vasco16.(param), 'min_max', min_max);
%     
%     % cut values above/below min/max limits
%     x = mio_min_max_cut(x, opt.vasco16.(param).min_max);
    
    % write file
    mdm_nii_write(x, fn{c}, h);
end


