function opt = dki_pa_opt(opt)
% function opt = dki_pa_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dki_pa.present = 1;

opt.dki_pa = msf_ensure_field(opt.dki_pa, ...
    'do_heteroscedasticity_correction', 1);

opt.dki_pa = msf_ensure_field(opt.dki_pa, ...
    'do_include_b_tensor_anisotropy', 0);

% control which maps that are generated as nifti files
opt.dki_pa = msf_ensure_field(opt.dki_pa, 'fig_maps', {'s0','MD', 'MK'});

opt.dki_pa = msf_ensure_field(opt.dki_pa, 'fig_prefix', 'dki_pa');

opt.dki_pa = msf_ensure_field(opt.dki_pa, 'pa_method', 1);

% produce standard fa-color map
opt.dki_pa = msf_ensure_field(opt.dki_pa, 'fig_cmaps',{});
opt.dki_pa = msf_ensure_field(opt.dki_pa, 'fig_ccol', {});
opt.dki_pa = msf_ensure_field(opt.dki_pa, 'fig_ccolnorm',{});


