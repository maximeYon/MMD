function M = mio_mask_thresh(I, opt)
% function M = mio_mask_thresh(I, opt)

opt.mask.present = 1;
opt.mask = msf_ensure_field(opt.mask, 'b0_ind', 1);
opt.mask = msf_ensure_field(opt.mask, 'thresh', 0.1);


% define the mask from threshold of b0 image
I0 = I(:,:,:,opt.mask.b0_ind);
M = double(I0)/double(max(I0(:))) > opt.mask.thresh;

M = mio_mask_fill(M);
M = mio_mask_keep_largest(M);

