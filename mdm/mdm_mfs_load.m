function mfs = mdm_mfs_load(mfs_fn)
% function mfs = mdm_mfs_load(mfs_fn)

load(mfs_fn,'mfs');

if (~exist('mfs','var'))
    error('mfs not found after load');
end