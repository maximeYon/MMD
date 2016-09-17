function m = ndi16_1d_data2fit(signal, xps, opt, ind)
% function m = ndi16_1d_data2fit(signal, xps, opt, ind)

if (nargin < 3), opt = ndi16_opt(); end
if (nargin < 4), ind = ones(size(signal)) > 0; end

unit_to_SI = [max(signal) 1 1 1e-9 1e-9];

    function m = t2m(t) % convert local params to outside format
        
        % define model parameters
        s0          = t(1);
        v_int       = t(2);
        v_csf       = t(3);
        lambda      = 1.7;
        md_csf      = 3.0;

        m = [s0 v_int v_csf lambda md_csf] .* unit_to_SI;
    end
                            
    function s = my_1d_fit2data(t,varargin)
        
        m = t2m(t);

        % signal
        s = ndi16_1d_fit2data(m, xps);
        s = s(ind);
                
    end

% first fit constrained to oblates -0.5<d_delta<0
% guesses and bounds
t_lb      = [eps    eps  eps  ];
t_ub      = [2-eps 1-eps 1-eps];

ss = inf;
for c = 1:opt.ndi16.n_rep
    t_guess   = t_lb + rand(size(t_lb)) .* (t_ub - t_lb);
    
    t_tmp = lsqcurvefit(@my_1d_fit2data, t_guess, [], signal(ind), t_lb, t_ub, opt.ndi16.lsq_opts);
    
    s_tmp = my_1d_fit2data(t_tmp);
    
    ss_tmp = sum(  (s_tmp(ind) - signal(ind)).^2 );
    if (ss_tmp < ss)
        ss = ss_tmp;
        t = t_tmp;
    end
end
m = t2m(t);

end