
function m = dki_pa_1d_data2fit(signal, xps, opt, ind)
% function m = dki_pa_1d_data2fit(signal, xps, opt, ind)
%
% fit the DKI model to powder averaged data (geo or ari)

if (nargin < 4), ind = ones(size(signal)) > 0; end

% Exclude data points with zero or negative values
ind = ind & (signal > 0);

% Setup regressors 
b0 = ones(size(xps.b(ind),1), 1);
b2 = xps.b(ind,:)      * 1e-9 ;   % SI unit conversion
b4 = xps.b(ind,:).^2   * 1e-18;   % SI unit conversion

X = [b0 -b2 1/2 * b4];
    
if (opt.dki_pa.do_heteroscedasticity_correction)
    C2 = diag(signal(ind));
else
    C2 = 1;
end

if (isfield(xps, 'pa_w'))
    C2 = C2 .* diag(xps.pa_w);
end

tmp = (X' * C2 * X);

if (rcond(tmp) > 1e-10) % some small number
 
    % perform regression to estimate model parameters m
    m = tmp \ X' * C2 * real(log(signal(ind)));
    
    m(1) = exp(m(1));
    m(2) = m(2) * 1e-9;   %converting back to SI units
    m(3) = m(3) * 1e-18;  %converting back to SI units
    
else
    warning('rcond fail')
    m = zeros(1, size(X,2));
end