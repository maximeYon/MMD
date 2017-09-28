function fn = msf_find_fn(bp, pattern)
% function fn = msf_find_fn(bp, pattern)

d = dir(fullfile(bp, pattern));

if (numel(d) == 0)
    error('no match');
    
end

if (numel(d) ~= 1)
    for c = 1:numel(d)
        disp(d(c).name);
    end
        
    error('non-unique pattern');
end

fn = fullfile(bp, d(1).name);