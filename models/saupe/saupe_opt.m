function opt = saupe_opt(opt)
% function opt = saupe_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.saupe.present = 1;

opt.saupe = msf_ensure_field(opt.saupe, 'tmp', 1); 
opt.saupe = msf_ensure_field(opt.saupe, 'fig_maps', {'op'});
opt.saupe = msf_ensure_field(opt.saupe, 'fig_prefix', 'saupe');
opt.saupe = msf_ensure_field(opt.saupe, 'fig_cmaps',{'mask'});
opt.saupe = msf_ensure_field(opt.saupe, 'fig_ccol',{'t1x6prim'});
opt.saupe = msf_ensure_field(opt.saupe, 'fig_ccolnorm',{'mask'});
