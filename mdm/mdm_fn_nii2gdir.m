function gdir_fn = mdm_fn_nii2gdir(nii_fn)
% function gdir_fn = mdm_fn_nii2gdir(nii_fn)
%
% convert a filename ending in .nii or .nii.gz to one ending in _gdir.txt

[a,b,~] = fileparts(nii_fn);

% handle possible .nii.gz extension
[~,b,~] = fileparts(b); 

gdir_fn = fullfile(a, [b '_gdir.txt']);

