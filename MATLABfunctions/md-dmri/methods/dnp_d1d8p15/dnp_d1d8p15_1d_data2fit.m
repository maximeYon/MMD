function m = dnp_d1d8p15_1d_data2fit(signal, xps, opt, ind)
% function m = dnp_d1d8p15_1d_data2fit(signal, xps, opt, ind)

if (nargin < 4), ind = (signal > 0.1 * max(signal)); end

unit_to_SI = [max(signal) 1 1 1e3 1e4 1 1e2];

    function m = t2m(t) % convert local params to outside format
        
        % define model parameters
        s0 = t(1);
        psd = t(2);
        rdnp = t(3);
        rsd = t(4);
        rcp = t(5);
        r1 = t(6);
        r1rho = t(7);

        m = [s0 psd rdnp rsd rcp r1 r1rho] .* unit_to_SI;
    end
                            
    function s = my_1d_fit2data(t,varargin)
        
        m = t2m(t);

        % signal
        s = dnp_d1d8p15_1d_fit2data(m, xps);
        s = s(ind);
                
    end

% first fit constrained to crosss -0.5<d_delta<0
% guesses and bounds
m_guess   = [max(signal)   .5 .2   5e3 2e4 .3   1e3];
m_lb      = [0             -1 0   0   0   0   0];
m_ub      = [10*max(signal) 1 1e1 1e6 1e6 1e1 1e4];

t_guess   = m_guess./unit_to_SI;
t_lb      = m_lb./unit_to_SI;
t_ub      = m_ub./unit_to_SI;

t = lsqcurvefit(@my_1d_fit2data, t_guess, [], signal(ind), t_lb, t_ub, opt.dnp_d1d8p15.lsq_opts);

m_cross = t2m(t);

signal_fit = dnp_d1d8p15_1d_fit2data(m_cross, xps);
chisq_cross = var(signal(ind) - signal_fit(ind));

% second fit constrained to diags 0<d_delta<1
% guesses and bounds
m_guess   = [max(signal)   -.5 .2   5e3 2e4 .3   1e3];
m_lb      = [0             -1 0   0   0   0   0];
m_ub      = [10*max(signal) 1 1e1 1e6 1e6 1e1 1e4];
                
t_guess   = m_guess./unit_to_SI;
t_lb      = m_lb./unit_to_SI;
t_ub      = m_ub./unit_to_SI;

t = lsqcurvefit(@my_1d_fit2data, t_guess, [], signal(ind), t_lb, t_ub,...
    opt.dnp_d1d8p15.lsq_opts);

m_diag = t2m(t);

signal_fit = dnp_d1d8p15_1d_fit2data(m_diag, xps);
chisq_diag = var(signal(ind) - signal_fit(ind));


if chisq_cross<chisq_diag
    m = m_cross;
else
    m = m_diag;
end


% if (opt.dnp_d1d8p15_.do_plot)
%     signal_fit = dnp_d1d8p15_1d_fit2data(m, xps);
%     semilogy(xps.b,signal,'.',xps.b,signal_fit,'o');
%     set(gca,'YLim',m(1)*[.01 1.2])
%     title(['D_\Delta = ' num2str(m(3),2)])
%     pause(0.05);
% end

end