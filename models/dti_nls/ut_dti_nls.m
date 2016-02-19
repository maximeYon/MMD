function fn = ut_dti_nls(c_ut)
% function fn = ut_dti_nls(c_ut)
%
% if c_ut is not supplied, the function returns the number of unit tests

% n_ut = number of unit tests
n_ut = 2;

if (nargin == 0)
    fn = n_ut;
    return;
end


mdiff = @(a,b,c) ~ ((a < (b + c)) && (a > (b - c)));

opt = dti_nls_opt;

xps.n = 13;
xps.bt = dtd_1x3_to_1x6(...
    [0 ones(1,12) * 1e9]', zeros(1,13)', [1 0 0; uvec_elstat_12dir]);


switch (c_ut)
    
    case 1 % make synthetic data, check fa and md
        fn = 'dti_nls_1d_data2fit.m';
        
        dt_3x3 = [2 0 0; 0 1 0; 0 0 0.3] * 1e-9;
        dt_1x6 = dtd_3x3_to_1x6(dt_3x3);
        
        m_exp = [1 dt_1x6];
        
        signal = exp(-dt_1x6 * xps.bt')';
        signal = signal + 0.0001 * randn(size(signal));
        
        m = dti_nls_1d_data2fit(signal, xps, opt);
        
        tol = [0.01 [1 1 1 1 1 1] * 1e-5];
        
        for c = 1:numel(m_exp)
            if (mdiff(m(c), m_exp(c), tol(c)))
                error('unexpected model param %i fit: %f/%f (r = %f)\n', c, m(c), m_exp(c), m(c) / (m_exp(c) + m(c))); 
            end
        end
        
    case 2
        fn = 'dti_nls_1d_fit2data.m';
        
        dt_3x3 = [2 0 0; 0 1 0; 0 0 0.3] * 1e-9;
        
        dt_1x6 = dtd_3x3_to_1x6(dt_3x3);
        
        m = [1 dt_1x6];
        s = dti_nls_1d_fit2data(m, xps);
        
        signal = exp(-dt_1x6 * xps.bt');
        
        for c = 1:numel(signal)
            if (mdiff(s(c), signal(c), 1e-6)), error('wrong model?'); end
        end
        
        
end