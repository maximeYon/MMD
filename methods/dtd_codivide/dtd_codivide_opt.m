function opt = dtd_codivide_opt(opt)
% function opt = dtd_codivide_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtd_codivide.present = 1;

opt.dtd_codivide = msf_ensure_field(opt.dtd_codivide, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e4));


opt.dtd_codivide = msf_ensure_field(opt.dtd_codivide, 'fig_maps', ...
    {'s0','v_at','v_fw','md_t'});

opt.dtd_codivide = msf_ensure_field(opt.dtd_codivide, 'fig_prefix', 'dtd_codivide');

opt.dtd_codivide = msf_ensure_field(opt.dtd_codivide, 'n_rep', 1);
