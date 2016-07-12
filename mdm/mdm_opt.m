function opt = mdm_opt(opt)
% function opt = mdm_opt(opt)
%
% Specifies default options 
1;


% opt.nii_ext - eiterh '.nii' or '.nii.gz' for compressed nii

% optional
% opt.i_range - control loop dimensions in 4d fit functions
% opt.j_range
% opt.k_range

if (nargin < 1), opt.present = 1; end

opt = msf_ensure_field(opt, 'nii_ext', '.nii.gz');
opt = msf_ensure_field(opt, 'do_overwrite', 1);
opt = msf_ensure_field(opt, 'verbose', 0);

opt = msf_ensure_field(opt, 'do_recon', 1);


opt = msf_ensure_field(opt, 'xps_merge_clear_s_ind', 0);

opt = msf_ensure_field(opt, 'xps_merge_rethrow_error', 1);
opt = msf_ensure_field(opt, 'pa_rethrow_error', 1);


opt = msf_ensure_field(opt, 'filter_sigma', 0);


opt = msf_ensure_field(opt, 'do_mask', 1);
opt = msf_ensure_field(opt, 'do_data2fit', 1);
opt = msf_ensure_field(opt, 'do_fit2param', 1);
opt = msf_ensure_field(opt, 'do_param2nii', 1);


opt = msf_ensure_field(opt, 'do_report_pdf', 1);
opt = msf_ensure_field(opt, 'do_maps_pdf', 1);


opt = msf_ensure_field(opt, 'do_pa_abs', 1);


opt.mdm.present = 1;
opt.mdm = msf_ensure_field(opt.mdm, 'mask_suffix', 'mask');
opt.mdm = msf_ensure_field(opt.mdm, 'txt_read_skip_comments', 0);

