function [t_1x6,w] = dtd_pa_dtd_to_t1x6(m, n)
% function t1x6 = dtd_pa_dta_to_t1x6(dtd, n)
%
% Convert internal representation to Voigt notated tensors pointing in 
% the 'n' direction
%

if (nargin < 2), n = [1 0 0]; end

[~,par,perp,w] = dtd_pa_dist2par(dtd_pa_m2dtd(m));

t_1x6 = zeros(numel(par), 6);
for c = 1:numel(par)
    t_1x6(c,:) = tm_1x3_to_1x6(par(c), perp(c), n);
end
    
