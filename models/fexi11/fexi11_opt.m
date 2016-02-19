function opt = fexi11_opt(opt)
% function opt = fexi11_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.fexi11.present = 1;

opt.fexi11 = msf_ensure_field(opt.fexi11, 'tmp', 1); 
opt.fexi11 = msf_ensure_field(opt.fexi11, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.fexi11 = msf_ensure_field(opt.fexi11, 'do_plot', 0);

opt.fexi11 = msf_ensure_field(opt.fexi11, 'mde_b2_ind_start', 1);