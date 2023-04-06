function M = mio_mask_erode_inplane(I, n)
% function M = mio_mask_erode_inplane(I)
%
% erode the mask by 1 voxel

if (nargin < 2)
    n = 3;
end

M = convn((double(I) > 0), ones(n,n), 'same') >= n^2;

