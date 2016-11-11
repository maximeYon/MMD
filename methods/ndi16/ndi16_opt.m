function opt = ndi16_opt(opt)
% function opt = ndi16_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.ndi16.present = 1;

opt.ndi16 = msf_ensure_field(opt.ndi16, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e4));

opt.ndi16 = msf_ensure_field(opt.ndi16, 'n_rep', 1); 
