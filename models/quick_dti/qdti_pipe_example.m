function qdti_pipe_example(s, o_path, opt)
% function qdti_pipe_example(nii_filename, o)
%
% s      - input structure
% o_path - output path

if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);

% Run the analysis
out_fn = fullfile(o_path, 'qdti.mat');
mfs_fn = qdti_4d_data2fit(s, out_fn, opt);

% Save parameter maps
fn = qdti_4d_fit2param(mfs_fn, o_path);

% Mask maps?



