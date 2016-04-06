function res = dtd_pa_pipe(s_pa, paths, opt)
% function res = dtd_pa_pipe(s_pa, paths, opt)
%

if (nargin < 3), opt.present = 1; end

opt = dtd_pa_opt(opt);

% Prepare: mask etc
if opt.do_mask
    opt.mask.b0_ind = 1;
    s_pa = mio_mask_thresh(s_pa, paths.mask_fn, opt);
else
    s_pa.mask_fn = paths.mask_fn;
end

% Run the analysis
if opt.do_fit
    res = dtd_pa_4d_data2fit(s_pa, paths.mfs.dtd_pa_primary_fn, opt);
end

% Convert from primary to derived model fit parameters
if opt.do_derivedparam
    res = dtd_pa_4d_fit2param(paths.mfs.dtd_pa_primary_fn, paths.mfs.dtd_pa_derived_fn);
end

% Save report pdf   
if opt.do_reportpdf
    res = dtd_pa_mkpdf(paths,opt);
end

% Save nifti parameter maps    
if opt.do_mapsnii
    res = mdm_mfs2nii(paths.mfs.dtd_pa_derived_fn, paths.maps, opt.dtd_pa, opt);
end

res = 1;
