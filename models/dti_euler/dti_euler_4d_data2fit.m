function mfs_fn = dti_euler_4d_data2fit(s, mfs_fn, opt)
% function mfs_fn = dti_euler_4d_data2fit(s, o, opt)

if (nargin < 3), opt = []; end

ind = opt.dti_euler.ind_start:s.xps.n;

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal) dti_euler_1d_data2fit(signal, xps, opt, ind);
mfs_fn = mio_fit_model(f, s, mfs_fn, opt);

