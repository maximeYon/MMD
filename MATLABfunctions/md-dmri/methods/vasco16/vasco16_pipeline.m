function fn = vasco16_pipeline(s, o_path, opt)
% function fn = vasco16_pipeline(s, o_path, opt)

% unit of v2 seems wrong
% calculation of alpha must be checked
warning('a few things must be checked, see this comment');


if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);

% First check that all fields are present, and init o_path
vasco16_check_xps(s.xps);
msf_mkdir(o_path);

% Prepare: mask, average and smooth
s = mio_mask_simple(s, o_path);
s = mio_powder_average(s, o_path, opt);
s = mio_smooth_4d(s, o_path, 0.5, opt);

% Run the analysis
out_fn = fullfile(o_path, 'vasco16.mat');

% limit the analysis to a square in one slice
if (1)
    opt.i_range = (30:80) + 10;
    opt.j_range = (30:80) + 10;
    opt.k_range = 12;
end

% Start the parallel pool
% parpool(4);

% Run the analysis
mfs_fn = vasco16_4d_data2fit(s, out_fn, opt);

% Save parameter maps
fn = vasco16_4d_fit2param(mfs_fn, o_path);
