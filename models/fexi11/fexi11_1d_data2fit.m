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

S0 = signal(xps.mde_b2_ind==1);

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



% exp fit - pre-fit
% tmp_guess = [signal(1) 1e-9];
% tmp_lb = [.9*signal(1) 1e-11];
% tmp_ub = [1.1*signal 3e-9];
% index = find(xps.mde_b1 == 0);
% Xin = xps.mde_b1(index);
% Xin = Xin(opt.nbstart:end);
% Yin = signal(index);
% Yin = Yin(opt.nbstart:end);

%Xin = xps.mde_b2;

%Xin = Xin(opt.nbstart:end);
% Yin = signal(index);
% Yin = Yin(opt.nbstart:end);

%%%%%%%%%%% for pre-fitting
% function Y = fexp(Pin,Xin,Pnorm,Xnorm,Ynorm);
%
% if nargin == 2
%     Pnorm = ones(size(Pin));
%     Xnorm = 1;
%     Ynorm = 1;
% end
%
% Pin = Pin.*Pnorm;
% Xin = Xin*Xnorm;
%
% Y0 = Pin(1);
% D = Pin(2);
%
% Y = Y0.*exp( - D.*Xin);
% Y = Y/Ynorm;
end



