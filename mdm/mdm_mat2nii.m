function mdm_mat2nii(I, path, sdim)
% function mdm_mat2nii(I, path, sdim)
%
% Makes nifti from I using path and header h

if ndims(I) == 2
    I = permute(I,[4 3 1 2]);
elseif ndims(I) == 3
    I = permute(I,[4 1 2 3]);
end

h = mdm_nii_h_empty;

h.pixdim(1+(1:length(sdim))) = sdim;

mdm_nii_write(I, path, h, 0)
    