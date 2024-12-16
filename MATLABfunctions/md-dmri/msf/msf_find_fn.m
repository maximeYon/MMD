function fn = msf_find_fn(bp, pattern, do_throw_error)
% function fn = msf_find_fn(bp, pattern)

if (nargin < 2), error('at least two input arguments required'); end
if (nargin < 3), do_throw_error = 1; end

% Allow returning empty filename if not found
if (do_throw_error == -1)
    f = @(x) 1; % be silent
    fn = [];
elseif (do_throw_error)
    f = @(x) error(x);
else
    f = @(x) disp(['msf_find_fn: ' x]);
    fn = []; % standard return
end

if (iscell(pattern))
    
    % Match to each pattern without throwing error if file not found
    fn = cell(1,numel(pattern));
    for c = 1:numel(pattern)
        fn{c} = msf_find_fn(bp, pattern{c}, -1);
    end
    
    % Throw error if not exactly one match
    ind = find(cellfun(@(x) ~isempty(x), fn));
    
    if (numel(ind) > 1) 
        f('Found multiple files, returning the first');
        fn = fn{ind(1)};
        return;
    elseif (numel(ind) == 0)
        f('msf_find_fn: Found no hits');
        fn = [];
        return;
    else
        fn = fn{ind};
        return;
    end
    
end

d = dir(fullfile(bp, pattern));

if (numel(d) == 0)
    f(sprintf('no match for pattern: %s', pattern)); return;
end

if (numel(d) ~= 1)
    str = []; for c = 1:numel(d), str = [str d(c).name ' ']; end
    f(sprintf('non-unique pattern (%s) with results: %s', pattern, str)); return;
end

fn = fullfile(bp, d(1).name);