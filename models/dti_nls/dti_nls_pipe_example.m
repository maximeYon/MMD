function fn = dti_nls_pipe_example(s, o_path, opt)
% function dti_nls_pipe_example(s, o_path, opt)
%
% s      - input structure
% o_path - output path

if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);


% Prepare: mask et c
s = mio_mask_simple(s, o_path);

% Run the analysis
out_fn = fullfile(o_path, 'dti_nls.mat');

opt.k_range = 30;

mfs_fn = dti_nls_4d_data2fit(s, out_fn, opt);

% Save parameter maps
fn = dti_nls_4d_fit2param(mfs_fn, o_path);
