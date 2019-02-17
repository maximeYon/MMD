function opt = dtd_opt(opt)
% function opt = dtd_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtd.present = 1;

opt.dtd = msf_ensure_field(opt.dtd, 'tmp', 1); 
opt.dtd = msf_ensure_field(opt.dtd, 'do_plot', 0);
opt.dtd = msf_ensure_field(opt.dtd, 'ind_start', 1);
opt.dtd = msf_ensure_field(opt.dtd, 'dmin', 1e-11);
opt.dtd = msf_ensure_field(opt.dtd, 'dmax', 5e-9);
opt.dtd = msf_ensure_field(opt.dtd, 'n_in', 1e2); % n_in: Number of nodes in NNLS inversion. [100 - 1000]
opt.dtd = msf_ensure_field(opt.dtd, 'n_out', 50);
opt.dtd = msf_ensure_field(opt.dtd, 'n_kill', 0);
opt.dtd = msf_ensure_field(opt.dtd, 'n_proliferation', 20);
opt.dtd = msf_ensure_field(opt.dtd, 'n_extinction', 20);
opt.dtd = msf_ensure_field(opt.dtd, 'ofuzz', .1*2*pi);
opt.dtd = msf_ensure_field(opt.dtd, 'dfuzz', .1);

opt.dtd = msf_ensure_field(opt.dtd, 'fig_maps', ...
    {'s0','mdiso','vdison','msqddelta','vsqddelta','cvdisosqddelta'});
opt.dtd = msf_ensure_field(opt.dtd, 'fig_prefix', 'dtd');
opt.dtd = msf_ensure_field(opt.dtd, 'fig_cmaps',{});
opt.dtd = msf_ensure_field(opt.dtd, 'fig_ccol',{});
opt.dtd = msf_ensure_field(opt.dtd, 'fig_ccolnorm',{});
opt.dtd = msf_ensure_field(opt.dtd, 'bin_disomax',[opt.dtd.dmax; 10^-9.5; opt.dtd.dmax; opt.dtd.dmax]);
opt.dtd = msf_ensure_field(opt.dtd, 'bin_disomin',[opt.dtd.dmin; opt.dtd.dmin; 10^-9.5; opt.dtd.dmin]);
opt.dtd = msf_ensure_field(opt.dtd, 'bin_dratiomax',[opt.dtd.dmax/opt.dtd.dmin; 10^.5; 10^.5; 10^-1]);
opt.dtd = msf_ensure_field(opt.dtd, 'bin_dratiomin',[10^.5; 10^-1; 10^-1; opt.dtd.dmin/opt.dtd.dmax]);

opt.dtd = msf_ensure_field(opt.dtd, 'odf_nnodes', 1000); %250, 500, 1000, 3994, or 15970
opt.dtd = msf_ensure_field(opt.dtd, 'odf_watsonkappa', 0.05);
