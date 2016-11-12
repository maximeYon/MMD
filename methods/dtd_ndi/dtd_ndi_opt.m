function opt = dtd_ndi_opt(opt)
% function opt = dtd_ndi_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtd_ndi.present = 1;

opt.dtd_ndi = msf_ensure_field(opt.dtd_ndi, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e4));

opt.dtd_ndi = msf_ensure_field(opt.dtd_ndi, 'n_rep', 1); 

opt.dtd_ndi = msf_ensure_field(opt.dtd_ndi, 'pa_method', 1);