function fn = ut_fexi11(c_ut)
% function fn = ut_fexi11(c_ut)
%
% if c_ut is not supplied, the function returns the number of unit tests

% n_ut = number of unit tests
n_ut = 1;

if (nargin == 0)
    fn = n_ut;
    return;
end

    function bool = approx_eq(a,b,c)
        bool = (a > (b - c)) & (a < (b + c));
    end


switch (c_ut)
    case 1
        
        xps.n            = 6;
        xps.mde_b1       = [0 0  1 1  1 1]' * 0.9e9;
        xps.mde_b2       = [0 1  0 1  0 1]' * 1.0e9;
        
        xps.mde_tm12     = [0 0  0 0  1 1]' * 0.3;
        
        xps.s_ind        = [1 1  2 2  3 3]';
        xps.mde_tm12_ind = [1 1  1 1  2 2]';
        xps.mde_b1_ind   = [1 1  2 2  2 2]';
        xps.mde_b2_ind   = [1 2  1 2  1 2]';
        
        mdm_xps_check(xps);
        
        adc = 0.7e-9;
        sigma = 0.5;
        axr = 4;
                
        m = [adc sigma axr 100 80 80];
        
        s = fexi11_1d_fit2data(m, xps);
        
        opt = fexi11_opt();
        
        m_fit = fexi11_1d_data2fit(s, xps, opt);
        
        tol = [1e-11 1e-3 1e-3 1e-3 1e-3 1e-3];
        for c = 1:numel(m)
            if (approx_eq(m_fit(c), m(c), tol(c)) == 0)
                error('something is strange in param %i', c); 
            end
        end
        
end

end