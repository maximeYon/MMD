function status = ut_dti_nls(c_ut)
% function status = ut_dti_nls(c_ut)
%
% if c_ut is not supplied, the function returns the number of unit tests

% n_ut = number of unit tests
n_ut = 1;

if (nargin == 0)
    status = n_ut;
    return;
end
    

switch (c_ut)
    case 1 % make synthetic phantom, check fa and md
        
        
end