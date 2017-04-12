function s = mdm_s_from_nii(nii_fn)
% function s = mdm_s_from_nii(nii_fn)
%
% converts a nii_fn to an input structure with two fields
%
% s.nii_fn
% s.xps
%
% assumes the xps filename can be constructred from the nii_fn

s.nii_fn = nii_fn;

if (exist(mdm_xps_fn_from_nii_fn(nii_fn), 'file'))
    s.xps = mdm_xps_load(mdm_xps_fn_from_nii_fn(nii_fn));
elseif (exist(mdm_fn_nii2gdir(nii_fn), 'file'))
    s.xps = mdm_xps_from_gdir(mdm_fn_nii2gdir(nii_fn));
end