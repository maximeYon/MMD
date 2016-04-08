function m = erf_1d_data2fit(signal, xps, opt, ind)
% function m = erf_1d_data2fit(signal, xps, opt, ind)

if (nargin < 3), opt = erf_opt(); end
if (nargin < 4), ind = ones(size(signal)) > 0; end

unit_to_SI = [max(signal) 1e-9 1];

    function m = t2m(t) % convert local params to outside format
        
        % define model parameters
        s0          = t(1);
        d_iso       = t(2);
        d_delta     = t(3);

        m = [s0 d_iso d_delta] .* unit_to_SI;
    end
                            
    function s = my_1d_fit2data(t,varargin)
        
        m = t2m(t);

        % signal
        s = erf_1d_fit2data(m, xps);
        s = s(ind);
                
    end

% first fit constrained to oblates -0.5<d_delta<0
% guesses and bounds
m_guess   = [max(signal)   1e-9  -0.25];
m_lb      = [0             1e-11 -0.49];
m_ub      = [2*max(signal) 3e-9  -0.001];

t_guess   = m_guess./unit_to_SI;
t_lb      = m_lb./unit_to_SI;
t_ub      = m_ub./unit_to_SI;

t = lsqcurvefit(@my_1d_fit2data, t_guess, [], signal(ind), t_lb, t_ub, opt.erf.lsq_opts);

m_oblate = t2m(t);

signal_fit = erf_1d_fit2data(m_oblate, xps);
chisq_oblate = var(signal(ind) - signal_fit(ind));

% second fit constrained to prolates 0<d_delta<1
% guesses and bounds
m_guess   = m_oblate.* [1 1 -1];
m_lb      = [0             1e-11 0];
m_ub      = [2*max(signal) 3e-9  1];
                
t_guess   = m_guess./unit_to_SI;
t_lb      = m_lb./unit_to_SI;
t_ub      = m_ub./unit_to_SI;

t = lsqcurvefit(@my_1d_fit2data, t_guess, [], signal(ind), t_lb, t_ub,...
    opt.erf.lsq_opts);

m_prolate = t2m(t);

signal_fit = erf_1d_fit2data(m_prolate, xps);
chisq_prolate = var(signal(ind) - signal_fit(ind));


if chisq_oblate<chisq_prolate
    m = m_oblate;
else
    m = m_prolate;
end


if (opt.erf.do_plot)
    signal_fit = erf_1d_fit2data(m, xps);
    semilogy(xps.b,signal,'.',xps.b,signal_fit,'o');
    set(gca,'YLim',m(1)*[.01 1.2])
    title(['D_\Delta = ' num2str(m(3),2)])
    pause(0.05);
end

end