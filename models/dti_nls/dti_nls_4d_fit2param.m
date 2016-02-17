function fn = dti_nls_4d_fit2param(mfs_fn, o_path, opt)
% function fn = dti_nls_4d_fit2param(mfs_fn, o_path, opt)
%
% Creates meaningful parameters from the model fit structure
%
% Input:
%
% mfs_fn - Path to .mat file with model fit structure
%
% o_path - Output path, where parameters maps are stored
%
% opt    - Options, optional argument
%
%
% Output:
%
% fn     - A cell array with paths to parameter maps that were written

if (nargin < 3), opt = []; end
    
opt = mdm_opt(opt);

mfs = mdm_mfs_load(mfs_fn);

    function o = my_core(m, do)
        
        if (~do)
            o = zeros(1,7);
        else
            
            C = [...
                m(2) m(5) m(7);
                0    m(3) m(6);
                0     0   m(4)];
            
            o = [m(1) dtd_3x3_to_1x6(C' * C)];
        end
    end

o = mio_volume_loop(@my_core, mfs.m, mfs.mask);

mfs.s0 = o(:,:,:,1); 
mfs.dt = o(:,:,:,2:7); 


% for now, use qdti's algorithm
fn = qdti_4d_fit2param(mfs, o_path, opt);


end