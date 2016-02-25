function fn = ut_uvec(c_ut)
% function fn = ut_uvec(c_ut)
%
% if c_ut is not supplied, the function returns the number of unit tests

% n_ut = number of unit tests
n_ut = 9;

if (nargin == 0)
    fn = n_ut;
    return;
end

switch (c_ut)
    
    case {1,2,3,4,5,6,7,8,9} 
        
        fn = {
            'uvec_elstat_6dir.m'
            'uvec_elstat_12dir.m'
            'uvec_elstat_32dir.m'
            'uvec_elstat_64dir.m'
            'uvec_elstat_128dir.m'
            'uvec_elstat_256dir.m'
            'uvec_dodeca.m'
            'uvec_icosa.m'
            'uvec_tricosa.m'};
        
        fn = fn{c_ut};
        
        u = feval(fn(1:end-2));
        
        if (size(u,2) ~= 3)
            error('%s: Second dimension not zero', fn);
        end
end