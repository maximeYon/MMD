function opt = vasco16_opt(opt)
% function opt = vasco16_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.vasco16.present = 1;

opt.vasco16 = msf_ensure_field(opt.vasco16, 'tmp', 1); 
opt.vasco16 = msf_ensure_field(opt.vasco16, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.vasco16 = msf_ensure_field(opt.vasco16, 'do_plot', 0);