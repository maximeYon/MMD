function mfs_fn = dtd_gamma_4d_data2fit(s, mfs_fn, opt)
% function mfs_fn = dtd_gamma_4d_data2fit(s, mfs_fn, opt)

if (nargin < 3), opt = []; end

ind = 1:s.xps.n;

% set b_eta = 0 if it is not in the xps already
if (~isfield(s.xps, 'b_eta')), s.xps.b_eta = zeros(size(s.xps.b)); end

% Verify the xps
dtd_gamma_check_xps(s.xps);

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal) dtd_gamma_1d_data2fit(signal, xps, opt, ind);

mfs_fn = mio_fit_model(f, s, mfs_fn, opt);

