function dps = codivide16_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = codivide16_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

    
opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them
dps.s0     = dps.m(:,:,:,1);
dps.v_at   = mio_min_max_cut( dps.m(:,:,:,2), [ 0 1 ]);
dps.v_fw   = mio_min_max_cut( dps.m(:,:,:,3), [ 0 1 ]);
dps.md_t   = mio_min_max_cut( dps.m(:,:,:,4) * 1e9, [0 4] );
dps.md_fw  = mio_min_max_cut( dps.m(:,:,:,5) * 1e9, [0 4] );


if (~isempty(dps_fn))
    mdm_dps_save(dps, dps.s, dps_fn, opt);
end

