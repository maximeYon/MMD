function mfs_fn = vasco16_4d_data2fit(s, o, opt)
% function mfs_fn = vasco16_4d_data2fit(s, o, opt)

if (nargin < 3), opt = []; end

% Verify the xps
vasco16_check_xps(s.xps);

% Loop over the volume and fit the model
opt         = vasco16_opt(opt);
xps         = s.xps;
ind_signal  = s.xps.b > 0;
f           = @(signal) vasco16_1d_data2fit(signal, xps, opt, ind_signal);
mfs_fn      = mio_fit_model(f, s, o, opt);

