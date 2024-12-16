function fn = ut_uvec(c_ut)
% function fn = ut_uvec(c_ut)
%
% if c_ut is not supplied, the function returns the number of unit tests

% n_ut = number of unit tests
n_ut = 4;

if (nargin == 0)
    fn = n_ut;
    return;
end

switch (c_ut)
    
    case {1,2,3}
        
        fn = {
            'uvec_dodeca.m'
            'uvec_icosa.m'
            'uvec_tricosa.m'};
        
        fn = fn{c_ut};
        
        u = feval(fn(1:end-2));
        
        if (size(u,2) ~= 3)
            error('%s: Second dimension not three', fn);
        end
        
    case 4
        fn = 'uvec_elstat.m';
        u = uvec_elstat(15);
        
        if (size(u,2) ~= 3)
            error('%s: Second dimension not three', fn);
        end
        
        
end