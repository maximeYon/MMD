function res = saupe_pipe(s, s_pa, paths, opt)
% function res = saupe_pipe(s, s_pa, paths, opt)
%

if (nargin < 4), opt.present = 1; end

% prepare options
opt = mdm_opt(opt);
opt = saupe_opt(opt);

% Prepare: mask etc
s = mio_mask_thresh(s, paths.mask_fn, opt);
s_pa.mask_fn = paths.mask_fn;
opt.do_mask = 0; % Use the same mask for dti, gamma, and erf

res = dti_euler_pipe(s, paths, opt);
res = gamma_pipe(s_pa, paths, opt);
res = erf_pipe(s_pa, paths, opt);

% Convert from primary to derived model fit parameters
if opt.do_derivedparam
    res = saupe_4d_fit2param(paths, paths.mfs.saupe_derived_fn);
end

% Save nifti parameter maps    
if opt.do_mapsnii
    res = mdm_mfs2nii(paths.mfs.saupe_derived_fn, paths.maps, opt.saupe, opt);
    % Save RGB parameter maps   
    res = mdm_mfs2cnii(paths.mfs.saupe_derived_fn, paths.maps, opt.saupe, opt);
end

res = 1;