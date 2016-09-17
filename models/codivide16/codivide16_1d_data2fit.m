function m = codivide16_1d_data2fit(signal, xps, opt, ind)
% function m = codivide16_1d_data2fit(signal, xps, opt, ind)
%
% s0     = m(1);
% v_at   = m(2);
% v_fw   = m(3);
% md     = m(4);
% md_fw  = m(5);
%

if (nargin < 3), opt = codivide16_opt(); end
if (nargin < 4), ind = ones(size(signal)) > 0; end

if (isfield(xps, 's_ind'))
    ns = max(xps.s_ind) - 1;
else
    ns = 0;
end

unit_to_SI = [max(signal) 1 1 1e-9 1e-9 ones(1,ns)];

    function m = t2m(t) % convert local params to outside format
        
        % define model parameters
        s0          = t(1);
        v_at        = t(2);
        v_fw        = t(3);
        md          = t(4);
        
        if (numel(t) >= (5 + ns)  )
            md_fw       = t(5);
        else
            md_fw   = 3.0;
        end
        
        sw          = t((end-(ns - 1)):end);
        
        m = [s0 v_at v_fw md md_fw sw] .* unit_to_SI;
    end

lambda = mean(signal);

    function s = my_1d_fit2data(t,varargin)
        
        m = t2m(t);
        s = codivide16_1d_fit2data(m, xps);
        
        penalty = exp(max(0, t(2) + t(3) - 1)) - 1;
        
        s = [s(ind); sum(penalty) * lambda];
        
    end

% first fit constrained to oblates -0.5<d_delta<0
% guesses and bounds
t_lb      = [0  eps  eps     eps 0.8 * ones(1,ns)];
t_ub      = [2 1-eps 1 - eps 1.8 1.2 * ones(1,ns)];

[t,ss_fit] = msf_fit(@my_1d_fit2data, [signal(ind); 0], t_lb, t_ub, opt.codivide16.n_rep, ...
    opt.codivide16.lsq_opts);

% test if a pure CSF fit is better
t_lb      = [0  1*eps  1.0-eps 1*eps 2.7 0.8 * ones(1,ns)];
t_ub      = [2  2*eps  1.0     2*eps 4.0 1.2 * ones(1,ns)];

[t2, ss_fit2] = msf_fit(@my_1d_fit2data, [signal(ind); 0], t_lb, t_ub, opt.codivide16.n_rep, ...
    opt.codivide16.lsq_opts);

if (ss_fit < ss_fit2)
    m = t2m(t);
else
    m = t2m(t2);
end


end