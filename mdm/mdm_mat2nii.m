function res = mdm_mat2nii(I,path,sdim)
%makes nifti from I using path and header h

if ndims(I) == 2
    I = permute(I,[4 3 1 2]);
end
if ndims(I) == 3
    I = permute(I,[4 1 2 3]);
end

h = mdm_nii_h_empty;

h.pixdim(1+(1:length(sdim))) = sdim;

mdm_nii_write(I, path, h, 0)
    
res = 1;
end

