function m = fexi11_1d_data2fit(signal, xps, opt, ind)
% function m = fexi11_1d_data2fit(signal, xps, opt, ind)
%
% Yields a 1xN vector 'm' with the fit parameters
% m(1)   - ADC0
% m(2)   - sigma
% m(3)   - AXR
% m(4:end) - vector of S0, length(m) = 3+max(s.xps.s_ind)

if (nargin < 4), ind = ~isnan(signal); end

signal = double(signal);

tmp_ind = xps.mde_b2_ind==min(xps.mde_b2_ind);
S0 = signal(tmp_ind);
S0(xps.s_ind(tmp_ind)) = S0;


t_guess    = [1.000  0.5 4.00 ones(size(S0'))];
t_lb       = [0.001  0.0 0.00 zeros(size(S0'))];
t_ub       = [4.000  1.0 20.0 2*ones(size(S0'))];

unit_to_SI = [1e-9   1.0 1.00  S0'];

    function s = fun(t,varargin)
        m = t .* unit_to_SI;
        s = fexi11_1d_fit2data(m, xps);
        s = s(ind);
    end

t = lsqcurvefit(@fun, t_guess, [], signal(ind), ...
    t_lb, t_ub, opt.fexi11.lsq_opts);

m = t .* unit_to_SI;


if (opt.fexi11.do_plot)
    signal_fit = fexi11_1d_fit2data(m, xps);
    x = (1:numel(signal))';
    plot(x,signal,'.',x,signal_fit,'o');
    pause(0.05);
end

end



