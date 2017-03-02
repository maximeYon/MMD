function opt = dki_lls_opt(opt)
% function opt = dki_lls_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dki_lls.present = 1;

opt.dki_lls = msf_ensure_field(opt.dki_lls, ...
    'do_heteroscedasticity_correction', 1);

% control which maps that are generated as nifti files
opt.dki_lls = msf_ensure_field(opt.dki_lls, 'fig_maps', ...
    {'s0','FA','MD', 'ad', 'rd', 'MK'});

opt.dki_lls = msf_ensure_field(opt.dki_lls, 'fig_prefix', 'dki_lls');

% produce standard fa-color map
opt.dki_lls = msf_ensure_field(opt.dki_lls, 'fig_cmaps',{'FA'});
opt.dki_lls = msf_ensure_field(opt.dki_lls, 'fig_ccol', {'u'});
opt.dki_lls = msf_ensure_field(opt.dki_lls, 'fig_ccolnorm',{'mask'});


