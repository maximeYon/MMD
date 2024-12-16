function M = mio_mask_cv(I, opt)
% function M = mio_mask_cv(I, opt)
%
% all images should have the same contrast
% 

opt.mask.present = 1;
opt.mask = msf_ensure_field(opt.mask, 'ind', 1:size(I,4));
opt.mask = msf_ensure_field(opt.mask, 'thresh', 0.2);


% define the mask from threshold of an image of the coefficient of
% variation
g = @(I) convn(I, ones(5,5,5), 'same');
f = @(I) g(std(I,[],4)) ./ g(mean(I,4));

CV = f(double(I(:,:,:,opt.mask.ind)));

M = CV < opt.mask.thresh;

M = mio_mask_fill(M);
M = mio_mask_keep_largest(M);

