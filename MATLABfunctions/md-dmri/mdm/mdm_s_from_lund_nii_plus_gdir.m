function s = mdm_s_from_lund_nii_plus_gdir(nii_fn)
% function s = mdm_s_from_lund_nii_plus_gdir(nii_fn)

s.nii_fn = nii_fn;
s.xps    = mdm_xps_from_gdir(mdm_fn_nii2gdir(nii_fn));
