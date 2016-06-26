function xps_fn = mdm_xps_fn_from_nii_fn(nii_fn)
% function xps_fn = mdm_xps_fn_from_nii_fn(nii_fn)

[path, name, ext] = msf_fileparts(nii_fn);

xps_fn = fullfile(path, [name '_xps.mat']);