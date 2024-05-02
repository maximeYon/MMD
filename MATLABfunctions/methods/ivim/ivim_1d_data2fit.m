function m = ivim_1d_data2fit(signal, xps, opt, ind)
% function m = ivim_1d_data2fit(signal, xps, opt, ind)
%
% s0      = m(1);
% f       = m(2);
% D_fast  = m(3);
% D_slow  = m(4);

if (nargin < 3), opt = ivim_opt(); end
if (nargin < 4), ind = ones(size(signal)) > 0; end

unit_to_SI = [max(signal) 1 1e-9 1e-9];

    function m = t2m(t) % convert local params to outside format
        
        % define model parameters
        s0          = t(1);
        f           = t(2);
        d_fast      = t(3); % 
        d_slow      = t(4);
        
        m = [s0 f d_fast d_slow] .* unit_to_SI;
    end

    function s = my_1d_fit2data(t,varargin)
        
        m = t2m(t);
        s = ivim_1d_fit2data(m, xps);
        
    end


% first fit: make it open
t_lb      = [0  eps  10 0.0];
t_ub      = [2 1-eps 100 3];

[t1,ss1] = msf_fit(@my_1d_fit2data, signal(ind), t_lb, t_ub, opt.ivim.n_rep, ...
    opt.ivim.lsq_opts);

% second fit: adc only
t_lb      = [0  eps  100 eps];
t_ub      = [2 1e-6  1000 3];

[t2,ss2] = msf_fit(@my_1d_fit2data, signal(ind), t_lb, t_ub, opt.ivim.n_rep, ...
    opt.ivim.lsq_opts);


if ( (ss2 / (numel(signal(ind)) - 2)) < (ss1 / (numel(signal) - 4)) )
    m = t2m(t2);
else
    m = t2m(t1);
end

% m = t2m(t2);
    



end