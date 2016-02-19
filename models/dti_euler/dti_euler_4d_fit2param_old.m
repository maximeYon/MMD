function fn = dti_euler_4d_fit2param(mfs_fn, o_path, opt)
% function fn = dti_euler_4d_fit2param(mfs_fn, o_path, opt)

if (nargin < 3), opt = []; end
    
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

% create parameter maps and save them
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


