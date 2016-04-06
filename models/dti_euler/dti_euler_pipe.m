function res = dti_euler_pipe(s, paths, opt)
% function res = dti_euler_pipe(s, paths, opt)
%

if (nargin < 3), opt.present = 1; end

opt = mdm_opt(opt);
opt = dti_euler_opt(opt);

% Prepare: mask etc
if opt.do_mask
    s = mio_mask_thresh(s, paths.mask_fn, opt);
else
    s.mask_fn = paths.mask_fn;
end

% Run the analysis
if opt.do_fit
    res = dti_euler_4d_data2fit(s, paths.mfs.dti_primary_fn, opt);
end

% Convert from primary to derived model fit parameters
if opt.do_derivedparam
    res = dti_euler_4d_fit2param(paths.mfs.dti_primary_fn, paths.mfs.dti_derived_fn);
end

% Save nifti parameter maps    
if opt.do_mapsnii
    res = mdm_mfs2nii(paths.mfs.dti_derived_fn, paths.maps, opt.dti_euler, opt);
    % Save RGB parameter maps   
    res = mdm_mfs2cnii(paths.mfs.dti_derived_fn, paths.maps, opt.dti_euler, opt);
end

res = 1;
