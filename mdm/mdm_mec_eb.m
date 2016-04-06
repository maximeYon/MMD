function s = mdm_mec_eb(s_target, s_source, p_fn, o_path, opt)
% function s = mdm_mec_eb(s_target, s_source, p_fn, o_path, opt)
% 
% Perform motion and eddy currect correction by registering to extrapolated
% references
%
% s_target  - input structure for target data
% s_source  - input structure for source data (should be low b-values)
% p_fn      - parameter filename, to elastix registration scheme
% o_path    - output path for the new files
% opt       - options (optional)
%
% Output
% s - s.nii_fn will be updated to refer to the corrected volume


% init
if (nargin < 4), opt.present = 1; end
opt = mdm_opt(opt);
msf_log(['Starting ' mfilename], opt);

[~,name]    = msf_fileparts(s_target.nii_fn);
ref_fn      = fullfile(o_path, [name '_ref' opt.nii_ext]);

% make reference
if (~exist(ref_fn, 'file') || opt.do_overwrite)
    msf_log('Making reference', opt);
    
    % load source data (do conventional MC first)
    [I,h] = mdm_nii_read(s_source.nii_fn);
    M = mdm_mask_load(s_source);
        
    % Create reference, improve fitting by some initial smoothing
    I     = mio_smooth_4d(single(I), 0.6);
    I_ref = mio_ref_extrapolate(I, s_source.xps, s_target.xps, M);
        
    % Rescale for easier comparisons
    I_true = mdm_nii_read(s_target.nii_fn);
    I_ref  = mio_rescale(I_ref, I_true, M);
        
    % Save output
    msf_delete(ref_fn);
    mdm_nii_write(I_ref, ref_fn, h);        
end

% Run coregistration
s = s_target;
s.nii_fn = mdm_coreg(s.nii_fn, ref_fn, p_fn, o_path, opt);

