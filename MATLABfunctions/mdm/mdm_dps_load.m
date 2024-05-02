function dps = mdm_dps_load(dps_fn)
% function dps = mdm_dps_load(dps_fn)
%
% load the derived parameter structure and check a few things

if (isstruct(dps_fn)) % assume it is a mfn struct already, error check below
    dps = dps_fn;
else
    load(dps_fn,'dps');
end

if (~exist('dps','var'))
    error('dps not found after load');
end

if (~isstruct(dps))
    error('dps must be a struct');
end

if (~isfield(dps, 'nii_h'))
    warning('dps.nii_h should exist and contain a nifti header');
end
