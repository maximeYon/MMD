function opt = dtd_covariance_opt(opt)
% function opt = dtd_covariance_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtd_covariance.present = 1;

opt.dtd_covariance = msf_ensure_field(opt.dtd_covariance, ...
    'do_heteroscedasticity_correction', 1);

% control which maps that are generated as nifti files
opt.dtd_covariance = msf_ensure_field(opt.dtd_covariance, 'fig_maps', ...
    {'s0','FA','C_MD', 'C_c', 'C_mu'});

opt.dtd_covariance = msf_ensure_field(opt.dtd_covariance, 'fig_prefix', 'dtd_covariance');

% produce standard fa-color map
opt.dtd_covariance = msf_ensure_field(opt.dtd_covariance, 'fig_cmaps',{'FA'});
opt.dtd_covariance = msf_ensure_field(opt.dtd_covariance, 'fig_ccol', {'u'});
opt.dtd_covariance = msf_ensure_field(opt.dtd_covariance, 'fig_ccolnorm',{'mask'});


