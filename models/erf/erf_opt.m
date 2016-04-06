function opt = erf_opt(opt)
% function opt = erf_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.erf.present = 1;

opt.erf = msf_ensure_field(opt.erf, 'tmp', 1); 
opt.erf = msf_ensure_field(opt.erf, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e3));
opt.erf = msf_ensure_field(opt.erf, 'do_plot', 0);
opt.erf = msf_ensure_field(opt.erf, 'fig_maps', {'s0','iso','delta','logratio','mu2aniso','mvlambda','ufa','cmu'});
opt.erf = msf_ensure_field(opt.erf, 'fig_prefix', 'erf');
