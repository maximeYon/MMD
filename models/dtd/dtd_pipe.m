function res = dtd_pipe(s, paths, opt)
% function res = dtd_pipe(s_pa, paths, opt)
%

if (nargin < 3), opt.present = 1; end

opt = dtd_opt(opt);

% Prepare: mask etc
if opt.do_mask
    s = mio_mask_thresh(s, paths.mask_fn, opt);
else
    s.mask_fn = paths.mask_fn;
end

% Run the analysis
if opt.do_fit
    res = dtd_4d_data2fit(s, paths.mfs.dtd_primary_fn, opt);
end

% Convert from primary to derived model fit parameters
if opt.do_derivedparam
    res = dtd_4d_fit2param(paths.mfs.dtd_primary_fn, paths.mfs.dtd_derived_fn);
end

% Save report pdf   
if opt.do_reportpdf
    res = dtd_mkpdf(paths,opt);
end

% Save nifti parameter maps    
if opt.do_mapsnii
    res = mdm_mfs2nii(paths.mfs.dtd_derived_fn, paths.maps, opt.dtd, opt);
end

res = 1;
