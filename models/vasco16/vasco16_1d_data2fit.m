function m = vasco16_1d_data2fit(signal, xps, opt, ind)
% function m = vasco16_1d_data2fit(signal, xps, opt, ind)

if (nargin < 4), ind = ones(size(signal)); end

unit_to_SI = [max(signal) 1e-6 1 1e-9 1e-9 1];

    function m = t2m(t) % convert local params to outside format
        
        % define model parameters
        s0          = t(1);
        D_tissue    = t(2); 
        vp          = t(3);
        f_blood     = t(4);
        D_blood     = 1.75; % AA 2016
        w           = t(5);
        
        m = [s0 vp f_blood D_blood D_tissue w] .* unit_to_SI;
    end

    function s = my_1d_fit2data_with_penalty(t,varargin)
        
        m = t2m(t);

        % calculate smooth heaviside regularization term
        vp = m(2) / unit_to_SI(2);
        nterm       = xps.n/30;
        lambda1     = 0.5;
        sigma1_neg  = 0.75;
        sigma1_pos  = 2.75;
        delta       = 0.4; % smoothness of transition
        gaussterm1  = 1-0.5*(tanh((vp-sigma1_neg)/delta) - ...
            tanh((vp-sigma1_pos)/delta));
        regterm1    = lambda1.*gaussterm1.*nterm;

        % make a signal plus a reg term
        s = vasco16_1d_fit2data(m, xps);
        
        s = [s(ind); regterm1 * mean(signal)];
        
    end

% S0, D, vd, f
t_guess   = [1  1    2   0.05 1.0];
t_lb      = [0  0    0   0.00 0.5];
t_ub      = [2 10  100   1.00 1.5];

% perform the fit
t = lsqcurvefit(@my_1d_fit2data_with_penalty, t_guess, [], [signal(ind);0], ...
    t_lb, t_ub, opt.vasco16.lsq_opts);

m = t2m(t);

end