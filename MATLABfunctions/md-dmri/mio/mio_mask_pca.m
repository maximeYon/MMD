function M = mio_mask_pca(I, opt)
% function M = mio_mask_pca(I, opt)
%
% pca-based filtering, inspired by an MRM paper (to do: find ref)

if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);
opt = mio_opt(opt);

sz = [size(I,1) size(I,2) size(I,3) size(I,4)];

I = reshape(single(I), prod(sz(1:3)), sz(4));
I = I - repmat(mean(I,2), 1, size(I,2));

[a,b] = eig(double(I' * I));

ind = (linspace(0,1,size(a,2))) > min(1 - 6 / size(a,2), 0.85);

I2 = I * a(:,ind) * a(:,ind)';

I2 = reshape(sqrt(mean( (I - I2).^2,2)), sz(1:3)) / mean(~ind);

I2 = imfilter(single(I2), ones(3,3,3)/27);

M = reshape(sqrt(mean(I.^2,2)),sz(1:3));
M = imfilter(M, ones(3,3,3) / 27);

M = (M ./ I2) > opt.mask.pca_threshold;


M = imfilter(single(M), ones(3,3,3)) > 10; % get rid of lonely voxels
M = imfilter(single(M), ones(3,3,3)) > 10; % get rid of lonely voxels

M = mio_mask_fill(M);
M = mio_mask_keep_largest(M);

