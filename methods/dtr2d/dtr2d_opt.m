function opt = dtr2d_opt(opt)
% function opt = dtr2d_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtr2d.present = 1;

opt.dtr2d = msf_ensure_field(opt.dtr2d, 'tmp', 1); 
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'do_plot', 0);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'ind_start', 1);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'dmin', 1e-11);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'dmax', 5e-9);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'r2min', 1);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'r2max', 20);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_in', 2e2); % n_in: Number of nodes in NNLS inversion. [100 - 1000]
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_out', 10);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_kill', 1);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_proliferation', 20);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_extinction', 20);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'ofuzz', .1*2*pi);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'dfuzz', .1);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'r2fuzz', .1);

opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_maps', ...
    {'s0','f_srfd','f_frfd','f_frsd','f_srsd','f_mrsds','f_mrsdp', 'f_mrsdl','f_residue',...
    'miso','viso_n','msaniso_n','vsaniso_n','mr2','vr2_n','cvisosaniso_n','cvisor2_n','cvsanisor2_n','ufa','fa','op'});
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_prefix', 'dtr2d');
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_cmaps',{'fa','cl','cp','ufa','ucl','ucp'});
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_ccol',{'t1x6','lambda33vec','lambda11vec','s1x6prim','s1x6prim','s1x6prim'});
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_ccolnorm',{'lambda33','mask','mask','slambda33prim','slambda33prim','slambda33prim'});

