function opt = gamma_opt(opt)
% function opt = gamma_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.gamma.present = 1;

opt.gamma = msf_ensure_field(opt.gamma, 'tmp', 1); 
opt.gamma = msf_ensure_field(opt.gamma, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e3));
opt.gamma = msf_ensure_field(opt.gamma, 'do_plot', 0);
opt.gamma = msf_ensure_field(opt.gamma, 'do_weight', 1);
opt.gamma = msf_ensure_field(opt.gamma, 'weight_sthresh', .2);
opt.gamma = msf_ensure_field(opt.gamma, 'weight_wthresh', 5);
opt.gamma = msf_ensure_field(opt.gamma, 'weight_mdthresh', 1e-9);
opt.gamma = msf_ensure_field(opt.gamma, 'fig_maps', ...
    {'s0','iso','mu2iso','mu2aniso','mvlambda','ufa','ciso','cmu'});
opt.gamma = msf_ensure_field(opt.gamma, 'fig_prefix', 'gamma');
