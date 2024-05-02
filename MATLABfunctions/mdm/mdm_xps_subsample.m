function xps = mdm_xps_subsample(xps, ind, opt)
% function xps = mdm_xps_subsample(xps, ind)

if (nargin < 3), opt = []; end

opt = mdm_opt(opt);

% help the user
if (opt.mdm_xps_throw_error)
    mdm_xps_check(xps);
end

if (numel(ind) ~= xps.n)
    error('expected logical array with %i elements ', xps.n);
end

f = fieldnames(xps);

for c = 1:numel(f)
    
    switch (f{c})
        
        case {'intent', 'n', 'c_volume'}
            1;
            
        case 's_ind'
            
            % Having non-continuous s_index is bad for fitting that allows
            % variable baseline signal for different series. Here we keep
            % the order of series, but make sure that no gaps appear.
            s_ind = xps.(f{c})(ind,:);
            us    = unique(s_ind)';
            u_ind = bsxfun(@eq, s_ind, us);
            s_new = mtimes(1:numel(us), u_ind')';
            
            xps.(f{c}) = s_new;
            
        otherwise
            
            try
                xps.(f{c}) = xps.(f{c})(ind,:);
            catch me
                if (opt.mdm_xps_throw_error)
                    disp(f{c});
                    rethrow(me);
                end
            end
            
    end
    
end

if (isfield(xps, 'b'))
    xps.n = numel(xps.b);
elseif (isfield(xps, 'u'))
    xps.n = size(xps.u, 1);
end
