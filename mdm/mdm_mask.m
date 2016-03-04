function s = mdm_mask(s, o, mask_fun, opt)
% function s = mdm_mask(s, o, mask_fun, opt)
%
% Mask the data in s.nii_fn using mask_fun (e.g. mio_mask_simple)
%
% Save the mask as a nifti to s.mask_fn. This field is created if it does
% not exist already

% construct the filename
if (~isfield(s.mask_fn))
    [~,name] = msf_fileparts(s.nii_fn);
    s.mask_fn = fullfile(o, [name opt.mdm.mask_suffix opt.nii_ext]);
end

if (exist(s.mask_fn, 'file') && (~opt.do_overwrite))
    return; 
end

% write the mask, don't care if we overwrite anything
[I,h] = mdm_nii_read(s.nii_fn);
M = mask_fun(I, opt);
mdm_nii_write(uint8(M), s.mask_fn, h);