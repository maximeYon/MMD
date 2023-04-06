function opt = dti_nls_opt(opt)
% function opt = dti_nls_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dti_nls.present = 1;

opt.dti_nls = msf_ensure_field(opt.dti_nls, 'tmp', 1); 
opt.dti_nls = msf_ensure_field(opt.dti_nls, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.dti_nls = msf_ensure_field(opt.dti_nls, 'do_plot', 0);