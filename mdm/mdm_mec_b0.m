function s = mdm_mec_b0(s, p_fn, o_path, opt)
% function s = mdm_mec_b0(s, p_fn, o_path, opt)
% 
% Perform motion and eddy currect correction by registering to b0
%
% s      - input structure
% p_fn   - parameter filename, to elastix registration scheme
% o_path - output path for the new files
% opt    - options (optional)
%
% Output
% s - s.nii_fn will be updated to refer to the corrected volume

if (nargin < 3), o_path = fileparts(s.nii_fn); end
if (nargin < 4), opt.present = 1; end
opt = mio_opt(opt);
msf_log(['Starting ' mfilename], opt);

% build file names
[~,name] = msf_fileparts(s.nii_fn);
ref_fn = fullfile(o_path, [name '_ref' opt.nii_ext]);

% save the volume with lowest b0 as the reference
if (~exist(ref_fn, 'file') || (opt.do_overwrite))
    [I,h] = mdm_nii_read(s.nii_fn);
    [~,c_b0] = min(s.xps.b);
    mdm_nii_write(I(:,:,:,c_b0), ref_fn, h);
end


s.nii_fn = mdm_coreg(s.nii_fn, ref_fn, p_fn, o_path, opt);
