function s = mio_mec_b0(s, p_fn, o_path, opt)
% function s = mio_mec_b0(s, p_fn, o_path, opt)
% 
% Perform motion and eddy currect correction by registering to b0

opt = mio_opt(opt);

% build file names
[~,name] = msf_fileparts(s.nii_fn);
ref_fn = fullfile(o_path, [name '_ref' opt.nii_ext]);


[I,h] = mdm_nii_read(s.nii_fn);


% for now, assume b0 is the first
c_b0 = 1;
mdm_nii_write(I(:,:,:,c_b0), ref_fn, h);


s.nii_fn = mio_coreg(s.nii_fn, ref_fn, p_fn, o_path, opt);


% here, we can do rotation of b-tensors et c according to transform
% parameters, the second argument of mio coreg