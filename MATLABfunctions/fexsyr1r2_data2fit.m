function m = fexsyr1r2_data2fit(signal, xps, opt, ind)
% function m = fexsyr1r2_data2fit(signal, xps, opt, ind)

if (nargin < 4), ind = ones(size(signal)) > 0; end

unit_to_SI = [max(signal(:))*[1 1] 1 1 100 100 1e-9*[1 1] 1];

    function m = t2m(t) % convert local params to outside format
        
        % define model parameters
        meqi          = t(1);
        meqe          = t(2);
        R1i           = t(3);
        R1e           = t(4);
        R2i           = t(5);
        R2e           = t(6);
        Di            = t(7);
        De            = t(8);
        kie           = t(9);
                
        m = [meqi meqe R1i R1e R2i R2e Di De kie] .* unit_to_SI;
    end

    function s = my_1d_fit2data(t,varargin)
        
        m = t2m(t);

        % signal
        s = fexsyr1r2_fit2data(m, xps);
        s = s(ind).*weight(ind);
                
    end

% weighting
weight = ones(xps.n,1);

% S0, D eigenvalues, Euler angles
m_guess   = [max(signal(:))*[1 1] 1 1 10 10 1e-11 1e-9 1];
m_lb      = [               0 0 0 0 0 0 0 0 0];
m_ub      = [100*max(signal(:))*[1 1] 10 10 500 500 5e-10 5e-9 500];
                
t_guess   = m_guess./unit_to_SI;
t_lb      = m_lb./unit_to_SI;
t_ub      = m_ub./unit_to_SI;

% perform the fit
t = lsqcurvefit(@my_1d_fit2data, t_guess, [], signal(ind).*weight(ind), t_lb, t_ub,...
    opt.fexsyr1r2.lsq_opts);

m = t2m(t);

% % optional second fit with weighting
% if opt.dti_euler.do_weight
%     signal_fit = dti_euler_1d_fit2data(m, xps);
% 
%     wthresh = opt.dti_euler.weight_wthresh;
%     sthresh = m(1)*opt.dti_euler.weight_sthresh;
%     %weight = .5*(erf(opt.dti_euler.weight_wthresh*(signal_fit/m(1) - opt.dti_euler.weight_sthresh)/opt.dti_euler.weight_sthresh)+1);
%     weight = .5*(erf(wthresh*(signal_fit - sthresh)/sthresh)+1);
% 
%     m_guess   = m;
%     t_guess   = m_guess./unit_to_SI;
% 
%     % perform the fit
%     t = lsqcurvefit(@my_1d_fit2data, t_guess, [], signal(ind).*weight(ind), t_lb, t_ub,...
%         opt.dti_euler.lsq_opts);
% 
%     m = t2m(t);
% end

% if (opt.dti_euler.do_plot)
%     signal_fit = dti_euler_1d_fit2data(m, xps);
%     [~,s_ind] = sort(signal_fit,'descend');
%     %semilogy(xps.b,signal,'.',xps.b,signal_fit,'o',xps.b,m(1)*weight,'x');
%     semilogy(1:xps.n,signal(s_ind),'.',1:xps.n,signal_fit(s_ind),'o',1:xps.n,m(1)*weight(s_ind),'x');
%     set(gca,'YLim',m(1)*[.01 1.2])
%     pause(0.05);
% end

end