function tpm = elastix_tpm_read(tpm_fn)
% function tpm = elastix_tpm_read(tpm_fn)
%
% Read a set of transform parameters from the Affine DTI Transform
%
% 

tpm = mdm_txt_read(tpm_fn);

tpm = cell2mat(cellfun(@(x) str2num(x)', tpm, 'uniformoutput', 0));

