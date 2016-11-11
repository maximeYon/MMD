function k = ut_mask_single_slice(mfs_fn, mask_fn)
% function k = ut_mask_single_slice(mfs_fn, mask_fn)

mfs = mdm_mfs_load(mfs_fn);

M = zeros(mfs.nii_h.dim(2:4)');
k = round(size(M,3)/2);
M(:,:,k) = 1;
M = M & mfs.mask;

mdm_nii_write(uint8(M), mask_fn, mfs.nii_h);
