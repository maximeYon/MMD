function msf_fn = ndi16_4d_data2fit(s, mfs_fn, opt)
% function msf_fn = ndi16_4d_data2fit(s, mfs_fn, opt)
%

if (nargin < 3), opt = []; end

% Loop over the volume and fit the model
xps     = s.xps; % this appears to improve parallel performance
f       = @(signal) ndi16_1d_data2fit(signal, xps, opt);
msf_fn  = mio_fit_model(f, s, mfs_fn, opt);

