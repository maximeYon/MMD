function s = mdm_mask(s, mask_fun, path, opt)
% function s = mdm_mask(s, mask_fun, path, opt)
%
% Mask the data in s.nii_fn using mask_fun (e.g. mio_mask_simple)
%
% Save the mask as a nifti to s.mask_fn. This field is created if it does
% not exist already

warning('legacy cody, use mdm_s_mask.m instead');

% init
if (nargin < 3) || (isempty(path)), path = fileparts(s.nii_fn); end
if (nargin < 4), opt.present = 1; end

s = mdm_s_mask(s, mask_fun, path, opt);

