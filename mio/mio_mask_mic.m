function M = mio_mask_mic(I, opt)
% function M = mio_mask_mic(I, opt)

if (nargin < 3), opt.present = 1; end
opt.mask.present = 1;
opt.mask = msf_ensure_field(opt.mask, 'b0_ind', 1);
opt.mask = msf_ensure_field(opt.mask, 'thres', 0.1);


% define the mask from 10% threshold of b0 image
I0 = I(:,:,:,opt.mask.b0_ind);
Imax = max(reshape(I0,[numel(I0) 1]));
M = I0/Imax > opt.mask.thresh;


M = mio_mask_fill(M,3);
M = mio_mask_fill(M,2);
M = mio_mask_fill(M,1);
M = mio_mask_fill(M,2);
M = mio_mask_fill(M,3);
M = mio_mask_keep_largest(M);
