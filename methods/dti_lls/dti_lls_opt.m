function opt = dti_lls_opt(opt)
% function opt = dti_lls_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dti_lls.present = 1;

opt.dti_lls = msf_ensure_field(opt.dti_lls, ...
    'do_heteroscedasticity_correction', 1);

% control which maps that are generated as nifti files
opt.dti_lls = msf_ensure_field(opt.dti_lls, 'fig_maps', ...
    {'s0','fa', 'md', 'ad', 'rd'});

opt.dti_lls = msf_ensure_field(opt.dti_lls, 'fig_prefix', 'dti_lls');

% produce standard fa-color map
opt.dti_lls = msf_ensure_field(opt.dti_lls, 'fig_cmaps',{'fa'});
opt.dti_lls = msf_ensure_field(opt.dti_lls, 'fig_ccol', {'u'});
opt.dti_lls = msf_ensure_field(opt.dti_lls, 'fig_ccolnorm',{'mask'});