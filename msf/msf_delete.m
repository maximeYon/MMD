function msf_delete(fn)

if (iscell(fn))
    for c = 1:numel(fn)
        msf_delete(fn{c});
    end
    return;
end

try
    if (isdir(fn) && exist(fn, 'dir'))
        rmdir(fn);
    elseif (exist(fn, 'file'))
        delete(fn);
    end
catch me
    disp(me.message);
end
