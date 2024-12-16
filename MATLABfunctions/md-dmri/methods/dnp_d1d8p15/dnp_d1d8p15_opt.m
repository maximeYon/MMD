function opt = dnp_d1d8p15_opt(opt)
% function opt = dnp_d1d8p15_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dnp_d1d8p15.present = 1;

opt.dnp_d1d8p15 = msf_ensure_field(opt.dnp_d1d8p15, 'tmp', 1); 
opt.dnp_d1d8p15 = msf_ensure_field(opt.dnp_d1d8p15, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e3));
opt.dnp_d1d8p15 = msf_ensure_field(opt.dnp_d1d8p15, 'do_plot', 0);
opt.dnp_d1d8p15 = msf_ensure_field(opt.dnp_d1d8p15, 'do_pa', 1);
opt.dnp_d1d8p15 = msf_ensure_field(opt.dnp_d1d8p15, 'fig_maps', {'s0','iso','saniso_n','ufa'});
opt.dnp_d1d8p15 = msf_ensure_field(opt.dnp_d1d8p15, 'fig_maps', {'s0','iso','delta','logratio','mu2aniso','vlambda','ufa','cmu'});
opt.dnp_d1d8p15 = msf_ensure_field(opt.dnp_d1d8p15, 'fig_prefix', 'dnp_d1d8p15');
