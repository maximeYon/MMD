function M = mio_mask_threshold(I, opt)
% function M = mio_mask_threshold(I, opt)

if (nargin < 2), opt.mask.present = 1; end
opt = mio_opt(opt);


% define the mask from threshold of b0 image
I0 = mean(I(:,:,:,opt.mask.b0_ind),4);
M = double(I0)/double(max(I0(:))) > opt.mask.threshold;

M = mio_mask_fill(M);
M = mio_mask_keep_largest(M);

