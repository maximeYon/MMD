function res = gamma_pipe(s_pa, paths, opt)
% function res = gamma_pipe(s_pa, paths, opt)
%

if (nargin < 3), opt.present = 1; end

opt = mdm_opt(opt);
opt = gamma_opt(opt);

% Prepare: mask etc
if opt.do_mask
    opt.mask.b0_ind = 1;
    s_pa = mio_mask_thresh(s_pa, paths.mask_fn, opt);
else
    s_pa.mask_fn = paths.mask_fn;
end

% Run the analysis
if opt.do_fit
    res = gamma_4d_data2fit(s_pa, paths.mfs.gamma_primary_fn, opt);
end

% Convert from primary to derived model fit parameters
if opt.do_derivedparam
    res = gamma_4d_fit2param(paths.mfs.gamma_primary_fn, paths.mfs.gamma_derived_fn);
end

% Save nifti parameter maps    
if opt.do_mapsnii
    res = mdm_mfs2nii(paths.mfs.gamma_derived_fn, paths.maps, opt.gamma, opt);
end

res = 1;
