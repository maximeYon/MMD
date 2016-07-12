function fn = ut_tensor_maths(c_ut)
% function fn = ut_tensor_maths(c_ut)
%
% if c_ut is not supplied, the function returns the number of unit tests

% 

% n_ut = number of unit tests
n_ut = 3;

if (nargin == 0)
    fn = n_ut;
    return;
end

switch (c_ut)
    
    case {1,2}
        fn = {'tm_1x6_to_3x3.m', 'tm_3x3_to_1x6.m'};
        fn = fn{c_ut};

        n = [ -0.4632   -0.5463    0.6978 ];
        t = tm_1x3_to_1x6(2,1,n);
        
        if any(abs(tm_3x3_to_1x6(tm_1x6_to_3x3(t)) - t) > eps)
            error('%s: Voigt back and forth not working', fn);
        end
   
    case 3
        fn = 'tm_1x6_eigvals.m';
        
        n = 100;
        t_1x6 = zeros(n, 6);
        L_supertrue = zeros(n, 3);
        for c = 1:n
            
            d = diag([0.1 1.0 2.0]);
            r = tm_euler_angles2rotmat(1*c,2*c,3*c);
            t_3x3 = r * d * r';
            t_1x6(c,:) = tm_3x3_to_1x6(t_3x3);
            
            L_supertrue(c,:) = sort(diag(d), 'descend');
        end
        
        % check quickly computed values versus those calculated with 'eigs'
        [L_true, U_true] = tm_1x6_eigvals(t_1x6, 2); % eigs
        [L_fast, U_fast] = tm_1x6_eigvals(t_1x6, 1); % Cardano
        
        for c = 1:n
            if any(abs( L_supertrue(c,:) - L_true(c,:)) > 1e-10)
                error('%s: True eigvals by  not working', fn);
            end
        end
        
        for c = 1:n
            if any(abs( L_fast(c,:) - L_true(c,:)) > 1e-10)
                error('%s: Quick eigvals by Cardano not working', fn);
            end
        end
        
        for c = 1:n
            if any( abs(abs(sum(U_fast(c,:) .* U_true(c,:))) - 1) > 1e-5 )
                error('%s: Quick eigvecs by Cardano not working', fn);
            end
        end
        
    case 4
        fn = 'tm_fa.m';
        
        dt_1x6 = [1 0 0 0 0 0; 0 0 0 0 0 0];
        
        fa = tm_fa(dt_1x6);
        
        if ( any( abs(fa - [1 0]') > eps ) )
            error('%s: FA calculation not working', fn);
        end
        
        
end