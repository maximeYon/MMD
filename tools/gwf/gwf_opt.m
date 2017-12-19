function opt = gwf_opt(opt)
% function opt = gwf_opt(opt)
%
% options structure for gwf analysis and plotting

opt.gwf.present = 1;

opt.gwf = msf_ensure_field(opt.gwf, ...
    'col', {[0.2 0.5 0.95], [0.4 0.9 0.1], [0.99 0.3 0.2]}); 


opt.gwf = msf_ensure_field(opt.gwf, 'nucleus', '1H');
opt.gwf = msf_ensure_field(opt.gwf, 'f_max', 100);