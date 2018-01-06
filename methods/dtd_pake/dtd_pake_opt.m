function opt = dtd_pake_opt(opt)
% function opt = dtd_pake_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtd_pake.present = 1;

opt.dtd_pake = msf_ensure_field(opt.dtd_pake, 'tmp', 1); 
opt.dtd_pake = msf_ensure_field(opt.dtd_pake, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e3));
opt.dtd_pake = msf_ensure_field(opt.dtd_pake, 'do_plot', 0);
<<<<<<< HEAD
opt.dtd_pake = msf_ensure_field(opt.dtd_pake, 'do_pa', 1);
opt.dtd_pake = msf_ensure_field(opt.dtd_pake, 'fig_maps', {'s0','iso','saniso_n','ufa'});
=======
opt.dtd_pake = msf_ensure_field(opt.dtd_pake, 'fig_maps', {'s0','iso','delta','logratio','mu2aniso','vlambda','ufa','cmu'});
>>>>>>> markus-nilsson/master
opt.dtd_pake = msf_ensure_field(opt.dtd_pake, 'fig_prefix', 'dtd_pake');
