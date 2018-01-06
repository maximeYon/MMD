function M = mio_mask_erode(I)
% function M = mio_mask_erode(I)
%
% erode the mask by 1 voxel

n = 3;
M = convn((double(I) > 0), ones(n,n,n), 'same') >= n^3;

