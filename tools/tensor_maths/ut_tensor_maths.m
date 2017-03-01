function fn = ut_tensor_maths(c_ut)
% function fn = ut_tensor_maths(c_ut)
%
% if c_ut is not supplied, the function returns the number of unit tests

% 

% n_ut = number of unit tests
n_ut = 5;

if (nargin == 0), fn = n_ut; return; end

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
        
        
    case 5
        fn = 'tm_tpars_to_1x6.m';
        
        % define b-value, direction (u) and b_delta
        b = 1.0;
        u = uvec_elstat(50);
        b_delta = linspace(-0.5, 1.0, size(u, 1))';

        % compute b-tensor from the tensor parameters above
        bt = tm_tpars_to_1x6(b, b_delta, u);

        % check that we get the same parameter back
        b_delta_out = zeros(size(b_delta));
        for c = 1:size(bt, 1)
            tp = tm_1x6_to_tpars(bt(c,:));
            b_delta_out(c) = tp.delta;
        end
        
        if (any( abs(b_delta - b_delta_out) > 1e-10))
            error('%s, ut_tensor_maths test %i, b_delta error', fn, c_ut);
        end            
        
        
end