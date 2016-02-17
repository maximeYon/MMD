function opt = mio_opt(opt)
% function opt = mio_opt(opt)

opt.mio.present = 1;

opt.mio = msf_ensure_field(opt.mio, 'n_core', 1);
