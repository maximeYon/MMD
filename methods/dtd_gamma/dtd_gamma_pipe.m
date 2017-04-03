function fn = dtd_gamma_pipe(s, paths, opt)
% function fn = dtd_gamma_pipe(s, paths, opt)
%
% s     - input structure
% paths - either a pathname or a path structure (see mdm_paths)
% opt   - (optional) options that drive the pipeline
%
% fn    - a cell arary with filenames to generated nii files
%
% Gamma fit
% Does powder averaging
% Lasic et al, Front. Phys. 2, 11 (2014).
% http://dx.doi.org/10.3389/fphy.2014.00011
% Modified for general b-tensor shapes


if (nargin < 3), opt.present = 1; end

opt = mdm_opt(opt);
opt = dtd_gamma_opt(opt);
paths = mdm_paths(paths, 'dtd_gamma');

msf_log(['Starting ' mfilename], opt);

% Prepare: mask and powder average
s = mdm_s_mask(s, @mio_mask_threshold, paths.nii_path, opt);

if (opt.dtd_gamma.do_pa)
    s = mdm_s_powder_average(s, paths.nii_path, opt);
end

% Run the analysis
mdm_data2fit(@dtd_gamma_4d_data2fit, s, paths.mfs_fn, opt);
mdm_fit2param(@dtd_gamma_4d_fit2param, paths.mfs_fn, paths.dps_fn, opt);

% Save nifti parameter maps
fn = mdm_param2nii(paths.dps_fn, paths.nii_path, opt.dtd_gamma, opt);
