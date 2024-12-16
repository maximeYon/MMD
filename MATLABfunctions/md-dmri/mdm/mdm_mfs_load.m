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
    path_nii = split(mfs_fn, filesep);
    path_nii = join(path_nii(1:end-4), filesep);
    path_nii = path_nii{1};
    path_nii = [path_nii filesep 'nii_xps' filesep 'data.nii.gz'];
    mfs.nii_h = mdm_nii_read_header(path_nii);
    if (~isfield(mfs, 'nii_h'))
        warning('mfs.nii_h should exist and contain a nifti header');
    end
end

if (~isfield(mfs, 'mask'))
    warning('mfs.mask should contain a mask');
end

if (~isfield(mfs, 'm'))
    warning('mfs.m is often the best variabel name for the model fit parameters');
end