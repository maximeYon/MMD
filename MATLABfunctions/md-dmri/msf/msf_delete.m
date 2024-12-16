function msf_delete(fn)
% function msf_delete(fn)
%
% Deletes the file fn, or multiple files if fn is a cell array

if (iscell(fn))
    for c = 1:numel(fn)
        msf_delete(fn{c});
    end
    return;
end

try
    if (isdir(fn) && exist(fn, 'dir'))
        rmdir(fn, 's');
    elseif (exist(fn, 'file'))
        delete(fn);
    end
catch me
    disp(me.message);
end
