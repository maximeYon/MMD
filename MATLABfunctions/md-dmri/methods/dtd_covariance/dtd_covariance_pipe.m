function fn = dtd_covariance_pipe(s, paths, opt)
% function fn = dtd_covariance_pipe(s, paths, opt)
%
% s     - input structure
% paths - either a pathname or a path structure (see mdm_paths)
% opt   - (optional) options that drive the pipeline
%         opt.mask.thresh = 0.1, you may want to adjust it 
%            
% fn    - a cell arary with filenames to generated nii files

if (nargin < 2), paths = fileparts(s.nii_fn); end
if (nargin < 3), opt.present = 1; end

% Init structures
opt   = mdm_opt(opt);
opt   = dtd_covariance_opt(opt);

paths = mdm_paths(paths, 'dtd_covariance');   

msf_log(['Starting ' mfilename], opt);

% Check that the xps is proper
dtd_covariance_check_xps(s.xps, opt);

% Prepare mask
if (opt.do_mask)
    s = mdm_s_mask(s, @mio_mask_threshold, [], opt);
end

% Smooth data
if (opt.filter_sigma > 0)
    s = mdm_s_smooth(s, opt.filter_sigma, fileparts(s.nii_fn), opt);
end

% Fit and derive parameters
if (opt.do_data2fit)
    mdm_data2fit(@dtd_covariance_4d_data2fit, s, paths.mfs_fn, opt);
end
if (opt.do_fit2param)
    mdm_fit2param(@dtd_covariance_4d_fit2param, paths.mfs_fn, paths.dps_fn, opt);
end

% Save nifti parameter maps    
if (opt.do_param2nii)
    fn = mdm_param2nii(paths.dps_fn, paths.nii_path, opt.dtd_covariance, opt); 
end

