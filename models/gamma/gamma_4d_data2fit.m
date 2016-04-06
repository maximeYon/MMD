function res = gamma_4d_data2fit(s, mfs_fn, opt)
% res = gamma_4d_data2fit(s, mfs_fn, opt)

if (nargin < 3), opt = []; end

res = -1;

ind = 1:s.xps.n;

%Verify the xps
%dti_euler_mic_check_xps(s.xps);

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal) gamma_1d_data2fit(signal, xps, opt, ind);
dummy = mio_fit_model(f, s, mfs_fn, opt);

res = 1;