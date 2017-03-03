function m = dki_lls_1d_data2fit(signal, xps, opt, ind)
% function m = dki_lls_1d_data2fit(signal, xps, opt, ind)
%
% fit the DKI model using linear least squares (as implemented in the
% dtd_covariance method)

if (nargin < 4), ind = ones(size(signal)) > 0; end


opt.dtd_covariance.do_heteroscedasticity_correction = ...
    opt.dki_lls.do_heteroscedasticity_correction;

opt.dtd_covariance.do_dki = 1;

m = dtd_covariance_1d_data2fit(signal, xps, opt, ind);

