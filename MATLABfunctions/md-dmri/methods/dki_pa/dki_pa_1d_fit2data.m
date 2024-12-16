function s = dki_pa_1d_fit2data(m, xps)
% function s = dki_pa_1d_fit2data(m, xps)
%
% m(1) - s0
% m(2) - md
% m(3) - v(d)

% define data from fitted parameters
if (numel(m) == 3)
    s = m(1) * exp(-xps.b * m(2) + 0.5 * xps.b.^2 * m(3));
elseif (numel(m) == 4)
    s = m(1) * exp(-xps.b * m(2) + 0.5 * xps.b.^2 * m(3) + 0.5 * xps.b.^2 .* xps.b_delta.^2 * m(4));    
else
    error('unknown');
end
