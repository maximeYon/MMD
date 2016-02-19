function s = msf_rmfield(s, f)
% function s = msf_rmfield(s, f)

if (iscell(f))
    for c = 1:numel(f)
       s = my_rmfield(s, f{c}); 
    end
    return;
end

if (isfield(s, f)), s = rmfield(s, f); end
