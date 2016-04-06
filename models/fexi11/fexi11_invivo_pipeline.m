function fn = fexi11_invivo_pipeline(s, o_path, opt)
% function fn = dti_nls_pipe_example(s, o_path, opt)
%
% s      - input structure
% o_path - output path

% init stuff
if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);
msf_mkdir(o_path);
fexi11_check_xps(s.xps);


% Motion and eddy current correction
n_iter  = 100;
p_fn    = elastix_p_3dof(fullfile(o_path, 'elastix_p.txt'), n_iter);
s       = mio_mec_b0(s, p_fn, o_path, opt);

% Mask and powder average
s = mio_powder_average(s, o_path, opt);

% limit the analysis to a square in one slice
if (opt.do_debug)
    opt.do_overwrite = 1;
end

s = mio_smooth_4d(s, o_path, 0.5, opt);
s = mio_mask_simple(s, o_path, opt);


% Run the analysis
out_fn  = fullfile(o_path, 'fexi11.mat');
mfs_fn  = fexi11_4d_data2fit(s, out_fn, opt);

% Save parameter maps
fn = fexi11_4d_fit2param(mfs_fn, o_path);
