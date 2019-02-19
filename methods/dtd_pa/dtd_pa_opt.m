function opt = dtd_pa_opt(opt)
% function opt = dtd_pa_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtd_pa.present = 1;

opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'tmp', 1); 
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'do_plot', 0);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'do_pa', 1);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'ind_start', 2);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'dmin', 1e-11);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'dmax', 5e-9);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'n_in', 1e2); % n_in: Number of nodes in NNLS inversion. [100 - 1000]
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'n_out', 50);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'n_kill', 0);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'n_proliferation', 20);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'n_extinction', 20);
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'dfuzz', .1);

opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'fig_maps', ...
    {'s0','miso','viso_n','msaniso_n','vsaniso_n'});
opt.dtd_pa = msf_ensure_field(opt.dtd_pa, 'fig_prefix', 'dtd_pa');
