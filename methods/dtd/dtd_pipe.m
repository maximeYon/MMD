function fn = dtd_pipe(s, paths, opt)
% function fn = dtd_pipe(s, paths, opt)
%
% s     - input structure
% paths - either a pathname or a path structure (see mdm_paths)
% opt   - (optional) options that drive the pipeline
%
% fn    - a cell array with filenames to generated nii files

fn = '';

if (nargin < 3), opt.present = 1; end

opt = mdm_opt(opt);
opt = dtd_opt(opt);
paths = mdm_paths(paths);
msf_log(['Starting ' mfilename], opt);

% Check that the xps is proper
dtd_check_xps(s.xps, opt);

% Prepare mask
if (opt.do_mask)
    s = mdm_s_mask(s, @mio_mask_threshold, [], opt);
end

% Smooth data
if (opt.filter_sigma > 0)
    s = mdm_s_smooth(s, opt.filter_sigma, fileparts(s.nii_fn), opt);
end

% Run the analysis
if (opt.do_data2fit)
    if opt.do_bootstrap
        msf_mkdir(fileparts(paths.ind_fn));
        ind = (opt.dtd.ind_start-0) + round(rand([s.xps.n-(opt.dtd.ind_start-1),1])*(s.xps.n-(opt.dtd.ind_start-0)));
        save(paths.ind_fn, 'ind');
%         ind_fn = mdm_ind_save(ind, paths.ind_fn);
%         load(paths.ind_fn)
        opt.bootstrap.ind = ind;
    end
    mdm_data2fit(@dtd_4d_data2fit, s, paths.mfs_fn, opt);
end

if (opt.do_datafit2chisq)
    chisq_fn = mio_datafit2chisq(@dtd_1d_fit2data, s, paths.mfs_fn, paths.chisq_fn, opt);
end

if (opt.do_fit2param)
    mdm_fit2param(@dtd_4d_fit2param, paths.mfs_fn, paths.dps_fn, opt);
end

% Save nifti parameter maps    
if (opt.do_param2nii)
    fn = mdm_param2nii(paths.dps_fn, paths.nii_path, opt.dtd, opt);
end





