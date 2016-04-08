function mfs_fn = dti_nls_4d_data2fit(s, o, opt)
% function mfs_fn = dti_nls_4d_data2fit(s, o, opt)
%
% Loops over a 4D volume to produce fit parameters with the DTI nls model
%
% Input: 
%
% s   - data structure
%          s.nii_fn - full path to nifti filename with data
%          s.xps    - experimental parameter structure
%
% o   - output folder   
% 
% opt - options structure, optional
%
% Output:
%
% mfs_fn - path to .mat file with the ModelFitStructure (mfs)

if (nargin < 3), opt = []; end

% Verify the xps
dti_nls_check_xps(s.xps);

% Loop over the volume and fit the model
opt     = dti_nls_opt(opt);
xps     = s.xps;
f       = @(signal) dti_nls_1d_data2fit(signal, xps, opt);
mfs_fn  = mio_fit_model(f, s, o, opt);

