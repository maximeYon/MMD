function fn = ut_fexi11(c_ut)
% function fn = ut_fexi11(c_ut)
%
% if c_ut is not supplied, the function returns the number of unit tests

% n_ut = number of unit tests
n_ut = 2;

if (nargin == 0)
    fn = n_ut;
    return;
end

    function bool = approx_eq(a,b,c)
        bool = (a > (b - c)) & (a < (b + c));
    end


xps.n            = 6;
xps.mde_b1       = [0 0  1 1  1 1]' * 0.9e9;
xps.mde_b2       = [0 1  0 1  0 1]' * 1.0e9;

xps.mde_tm12     = [0 0  0 0  1 1]' * 0.3;

xps.s_ind        = [1 1  2 2  3 3]';
xps.mde_tm12_ind = [1 1  1 1  2 2]';
xps.mde_b1_ind   = [1 1  2 2  2 2]';
xps.mde_b2_ind   = [1 2  1 2  1 2]';

mdm_xps_check(xps);

switch (c_ut)
    case 1
        fn = 'fexi11_1d_data2fit';
        
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
        
        
    case 2
        fn = 'fexi11_4d_data2fit';
        
        
        n_x = 4;
        I = zeros(n_x,n_x,3,xps.n);
        
        adc = 0.7e-9;
        sigma = 0.5;
        axr = 4;
                
        m = [adc sigma axr 100 80 80];
        
        for i = 1:size(I,1)
            for j = 1:size(I,2)
                for k = 1:size(I,3)
                    I(i,j,k,:) = fexi11_1d_fit2data(m, xps);
                end
            end
        end
        
        p = msf_tmp_path;
        
        s.nii_fn = fullfile(p, 'I.nii.gz');
        s.mask_fn = fullfile(p, 'M.nii.gz');
        s.xps = xps;
        
        h = mdm_nii_h_empty;
        
        mdm_nii_write(single(I), s.nii_fn, h);
        mdm_nii_write(ones(size(I,1), size(I,2), size(I,3)), s.mask_fn, h);
        
        maps_fn = fexi11_4d_data2fit(s, fullfile(p, 'fexi.mat'));
        
        rmdir(p, 's');
        
end

end