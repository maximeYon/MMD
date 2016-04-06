function fn = ut_mdm(c_ut)
% function fn = ut_mdm(c_ut)
%
% Run unit tests on the files in this package

if (nargin == 0), fn = 1; return; end

switch (c_ut)
    
    case 1
        fn = 'mdm_nii_write.m, mdm_nii_read.m';
        
        tmp_fn = fullfile(msf_tmp_path, 'a.nii');
        
        s = [1 16 16 4];
        
        a = randn(s) + 1i * randn(s);
        
        
        mdm_nii_write(a, tmp_fn);
        
        o = mdm_nii_read(tmp_fn);
        msf_delete(tmp_fn);
        
        if (...
                (size(a,1) ~= size(o,1)) || ...
                (size(a,2) ~= size(o,2)) || ...
                (size(a,3) ~= size(o,3)) || ...
                (size(a,4) ~= size(o,4)))
            
            error('%s, ut_mdm test %i, dimension error', fn, c_ut);
        end
        
        if (sum(abs(a(:) - o(:))) > eps)
            error('%s, ut_mdm test %i, read/write different', fn, c_ut);
        end
        
        
        
end
