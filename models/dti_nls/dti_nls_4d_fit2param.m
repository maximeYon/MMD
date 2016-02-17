function fn = dti_nls_4d_fit2param(mfs_fn, o_path, opt)
% function dti_4d_fit2param(mfs, s, o_path, opt)

if (nargin < 3) opt = []; end
    
opt = mdm_opt(opt);

mfs = mdm_mfs_load(mfs_fn);

sz = size(mfs.m);

mfs.s0 = zeros(sz(1:3));
mfs.dt = zeros([sz(1:3), 6]);

    function mfs = my_core(m, do)
        
        if (~do)
            m = zeros(1,7);
        end
        
        C = [...
            m(2) m(5) m(7);
            0    m(3) m(6);
            0     0   m(4)];
        
        mfs.s0 = m(1);
        mfs.dt = dtd_3x3_to_1x6(C' * C);
    end

mfs = mio_volume_loop(@my_core, mfs.m, mfs.mask, mfs);


% for now, use qdti's algorithm
fn = qdti_4d_fit2param(mfs, o_path, opt);


end