function res = dtor1r2d_4d_data2fit(s, mfs_fn, opt)
% function mfs_fn = dtor1r2d_4d_data2fit(s, o, opt)

if (nargin < 3), opt = []; end

res = -1;

if opt.do_bootstrap
    ind = opt.bootstrap.ind;
else
    ind = opt.dtor1r2d.ind_start:s.xps.n;
end

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal) dtor1r2d_1d_data2fit(signal, xps, opt, ind);
dummy = mio_fit_model(f, s, mfs_fn, opt);

res = 1;