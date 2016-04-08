function erf_4d_data2fit(s, mfs_fn, opt)
% function erf_4d_data2fit(s, mfs_fn, opt)
%
% Error function fit
% Assumes powder averaging and axisymmetric b-tensors
% Eriksson et al., J. Chem. Phys. 142, 104201 (2015).
% http://dx.doi.org/10.1063/1.4913502


if (nargin < 3), opt = []; end

ind = 1:s.xps.n;

%Verify the xps
%dti_euler_mic_check_xps(s.xps);

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal) erf_1d_data2fit(signal, xps, opt, ind);
dummy = mio_fit_model(f, s, mfs_fn, opt);

