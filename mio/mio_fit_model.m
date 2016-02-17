function mfs_fn = mio_fit_model(fun, s, o, opt)
% function mfs = mio_fit_model(s, fun, opt)

if (nargin < 3), error('all inputs required'); end

% Read and reformat data
[I,h]  = mdm_nii_read(s.nii_fn);
M      = mdm_mask_load(s, opt);

% Analyze and store output
mfs.m       = mio_volume_loop(fun, I, M);
mfs.mask    = M;
mfs.nii_h   = h;

% Save data
mfs_fn = mdm_mfs_save(mfs, s, o, opt);
