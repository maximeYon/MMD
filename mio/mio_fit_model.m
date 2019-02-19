function mfs_fn = mio_fit_model(fun, s, o_fn, opt)
% function mfs_fn = mio_fit_model(fun, s, o_fn, opt)
%
% Fits the model specified in 'fun' to the data referred to in 's', 
% and saves the output in 'o_fn', which should have '.mat' as its extension
%
% Typically, 'fun' is a local function defined within model_4d_data2fit

if (nargin < 3), error('first three inputs required'); end
if (nargin < 4), opt.present = 1; end

opt = mio_opt(opt);

% Read and reformat data
[I,h]  = mdm_nii_read(s.nii_fn);
M      = mdm_mask_load(s, opt);

h.scl_slope = 1;
h.scl_inter = 0;

% Disallow model fits to complex data
if (any(imag(I) ~= 0)), I = abs(I); end 

% Analyze and store output
mfs.m       = mio_volume_loop(fun, I, M, opt);
mfs.mask    = M;
mfs.nii_h   = h;

% Save data
mfs_fn = mdm_mfs_save(mfs, s, o_fn, opt);
