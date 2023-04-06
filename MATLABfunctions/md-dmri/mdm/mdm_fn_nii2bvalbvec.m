function [bval_fn, bvec_fn] = mdm_fn_nii2bvalbvec(nii_fn)
% function [bval_fn, bvec_fn] = mdm_fn_nii2bvalbvec(nii_fn)
%
% convert a filename ending in .nii or .nii.gz to two files ending in
% .bval and .bvec

[nii_path, nii_name] = msf_fileparts(nii_fn);
bval_fn = fullfile(nii_path, [nii_name '.bval']);
bvec_fn = fullfile(nii_path, [nii_name '.bvec']);
