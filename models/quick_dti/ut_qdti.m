function status = ut_qdti(c_ut)
% function ut_qdti(c_ut)

% n_ut = number of unit tests
n_ut = 1;

if (nargin == 0)
    status = n_ut;
    return;
end
    

switch (c_ut)
    case 1 % make synthetic phantom, check fa and md
        
end