function opt = mdm_opt(opt)

if (nargin < 1), opt.present = 1; end

opt = ensure_field(opt, 'nii_ext', '.nii.gz');