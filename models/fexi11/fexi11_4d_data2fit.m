function mfs_fn = fexi11_4d_data2fit(s, o, opt)
% function mfs_fn = fexi11_4d_data2fit(s, o, opt)
%
% Loops over a 4D volume to produce fit parameters with the FEXI nls model
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

opt = fexi11_opt(opt);

if (isfield(s.xps,'mde_b2_ind'))
    ind = s.xps.mde_b2_ind >= opt.fexi11.mde_b2_ind_start;
else
    ind = ~isnan(s.xps.b);
end

% Verify the xps
fexi11_check_xps(s.xps);

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal) fexi11_1d_data2fit(signal, xps, opt, ind);
mfs_fn = mio_fit_model(f, s, o, opt);

