function res = dtd_4d_data2fit(s, mfs_fn, opt)
% function mfs_fn = dtd_4d_data2fit(s, o, opt)

if (nargin < 3), opt = []; end

res = -1;

ind = 1:s.xps.n;

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal) dtd_pa_1d_data2fit(signal, xps, opt, ind);
dummy = mio_fit_model(f, s, mfs_fn, opt);

res = 1;