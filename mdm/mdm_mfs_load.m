function mfs = mdm_mfs_load(mfs_fn)
% function mfs = mdm_mfs_load(mfs_fn)
%
% Load model fit structure

if (isstruct(mfs_fn)) % assume it is a mfn struct already, error check below
    mfs = mfs_fn;
else
    load(mfs_fn,'mfs');
end

if (~exist('mfs','var'))
    error('mfs not found after load');
end

if (~isstruct(mfs))
    error('mfs must be a struct');
end

if (~isfield(mfs, 'nii_h'))
    warning('mfs.nii_h should exist and contain a nifti header');
end

if (~isfield(mfs, 'mask'))
    warning('mfs.mask should contain a mask');
end

if (~isfield(mfs, 'm'))
    warning('mfs.m is often the best variabel name for the model fit parameters');
end