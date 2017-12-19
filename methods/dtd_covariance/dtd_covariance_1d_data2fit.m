function [m,cond] = dtd_covariance_1d_data2fit(signal, xps, opt, ind)
% function m = dtd_covariance_1d_data2fit(signal, xps, opt, ind)
%
% Diffusion tensor distribution (DTD) modeling using cumulant expansion:
% This function calculates the first two terms of the cumulant expansion 
% from Q-space trajectory (QTI) data. 
%
% In the cumulant expansion
% - The first term corresponds to the mean of the tensor distribution, a second order
%   tensor, which is closely realted to the diffusion tensor in DTI
% - The second term is the covariance of the tensor distribution, a fourth order tensor
%
% The output from this function will be a model parameter vector of
% dimension 1 + 6 + 21 = 28
%
% The second output is the condition number of the matrix used in the
% inversion

if (nargin < 4), ind = ones(size(signal)) > 0; end

% log S = X * B (B --> m, our model parameters)
%
% Correct for heteroscedasticity
%
% C * log S = C * X * B
%
% inv(X' * C' * C * X) * X' * C' * C * log S = B
%
% Compute this whole step for each iteration. It is slow, but use the
% parallell toolbox if you are in a hurry


% Exclude data points with zero or negative values
ind = ind & (signal > 0);


% Setup regressors for diffusion tensor distribution (DTD) model
b0 = ones(xps.n, 1);
b0 = b0(ind);

b2 = xps.bt(ind,:)      * 1e-9 ;   %SI unit conversion

if (opt.dtd_covariance.do_dki)
    % the dtd_covariance model can be reduced to the DKI model by using
    % the directions (normalized, stored in xps.u) and the b-values 
    % to create a 1x15 4th order tensor which captures a subset of the
    % information available if using the full covariance model
    b4 = tm_1x3_to_1x15(xps.u(ind,:)) .* repmat(xps.b(ind).^2 * 1e-18, 1, 15);
else
    b4 = tm_1x6_to_1x21(b2); % SI unit conversion already done for b2
end

X = [b0 -b2 1/2 * b4];
    
if (opt.dtd_covariance.do_heteroscedasticity_correction)
    C2 = diag(signal(ind));
else
    C2 = 1;
end

% Compute condition number
tmp = (X' * C2 * X);
cond = rcond(tmp);

if (cond > 1e-10) % some small number
 
    % perform regression to estimate model parameters m
    m = tmp \ X' * C2 * real(log(signal(ind)));
    
    m(1)    = exp(m(1));
    m(2:7)  = m(2:7)  * 1e-9;   %converting back to SI units
    m(8:end) = m(8:end) * 1e-18;  %converting back to SI units
    
else
    warning('rcond fail in dtd_covariance_1d_data2fit')

    % number of paramters, S0 + 6 for mean + 21 for cov = 28
    m = zeros(1, size(X,2));
end

