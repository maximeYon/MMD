function fn = dtd_saupe_pipe(s, paths, opt)
% function fn = dtd_saupe_pipe(s, paths, opt)

if (nargin < 4), opt.present = 1; end

% prepare options
opt = mdm_opt(opt);
opt = dti_euler_opt(opt);
opt = gamma_opt(opt);
opt = erf_opt(opt);
opt = dtd_saupe_opt(opt);

% Setup paths
paths = mdm_paths(paths);
paths.mfs_dti_euler_fn = fullfile(fileparts(paths.mfs_fn), 'mfs_dti_euler.mat');
paths.mfs_gamma_fn     = fullfile(fileparts(paths.mfs_fn), 'mfs_gamma.mat');
paths.mfs_erf_fn       = fullfile(fileparts(paths.mfs_fn), 'mfs_erf.mat');

% Prepare: mask etc
s = mdm_mask(s, @mio_mask_threshold, [], opt);

% Run sub analyses
mdm_data2fit(@dti_euler_4d_data2fit, s, paths.mfs_dti_euler_fn, opt);
mdm_data2fit(@gamma_4d_data2fit,     s, paths.mfs_gamma_fn, opt);
mdm_data2fit(@erf_4d_data2fit,       s, paths.mfs_erf_fn, opt);


% Convert from primary to derived model fit parameters
mdm_fit2param(@dtd_saupe_4d_fit2param, ...
    {paths.mfs_dti_euler_fn, paths.mfs_gamma_fn, paths.mfs_erf_fn}, ...
    paths.dps_fn, opt);

% Save nifti parameter maps    
fn = mdm_param2nii(paths.dps_fn, paths.nii_path, opt.dtd_saupe, opt);

