function opt = codivide16_opt(opt)
% function opt = codivide16_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.codivide16.present = 1;

opt.codivide16 = msf_ensure_field(opt.codivide16, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e4));


opt.codivide16 = msf_ensure_field(opt.codivide16, 'fig_maps', ...
    {'s0','v_at','v_fw','md_t'});

opt.codivide16 = msf_ensure_field(opt.codivide16, 'fig_prefix', 'codivide16');

opt.codivide16 = msf_ensure_field(opt.codivide16, 'n_rep', 1);
