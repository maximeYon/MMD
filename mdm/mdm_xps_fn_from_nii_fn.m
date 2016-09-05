function xps_fn = mdm_xps_fn_from_nii_fn(nii_fn)
% function xps_fn = mdm_xps_fn_from_nii_fn(nii_fn)

if (iscell(nii_fn))
    xps_fn = cellfun(@mdm_xps_fn_from_nii_fn, nii_fn, 'uniformoutput', 0);
    return;
end

[path, name] = msf_fileparts(nii_fn);

xps_fn = fullfile(path, [name '_xps.mat']);