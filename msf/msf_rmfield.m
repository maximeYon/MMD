function s = msf_rmfield(s, f)
% function s = msf_rmfield(s, f)
%
% Removes field 'f' from 's' is possible. If 'f' is a cell array, all
% fields in 'f' are removed

if (iscell(f))
    for c = 1:numel(f)
       s = my_rmfield(s, f{c}); 
    end
    return;
end

if (isfield(s, f)), s = rmfield(s, f); end
