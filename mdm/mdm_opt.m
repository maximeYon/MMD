function opt = mdm_opt(opt)
% function opt = mdm_opt(opt) 


% opt.nii_ext - eiterh '.nii' or '.nii.gz' for compressed nii

% optional
% opt.i_range - control loop dimensions in 4d fit functions
% opt.j_range
% opt.k_range

if (nargin < 1), opt.present = 1; end

opt = msf_ensure_field(opt, 'nii_ext', '.nii.gz');
