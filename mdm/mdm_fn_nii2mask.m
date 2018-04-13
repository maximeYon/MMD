function mask_fn = mdm_fn_nii2mask(nii_fn, opt)
% function mask_fn = mdm_fn_nii2mask(nii_fn, opt)

if nargin < 2
    opt = mdm_opt();
end

[path,name] = msf_fileparts(nii_fn);
mask_fn = fullfile(path, [name '_' opt.mdm.mask_suffix opt.nii_ext]);