function out_path = mdm_bruker_spen2nii_xps(data_path, out_path, rps, opt)
% function out_path = mdm_bruker_spen2nii_xps(data_path, out_path, rps, opt)
%
%
% Works in progress

if (nargin < 5), opt.present = 1; end

opt = mdm_opt(opt);

opt = mdm_opt();

nii_temp_fn = fullfile(out_path, ['data_temp' opt.nii_ext]);
xps_temp_fn = fullfile(out_path, 'data_temp_xps.mat');

mdm_bruker_spen180dwwave_2dseq2nii(data_path, nii_temp_fn, rps);
mdm_bruker_spen180dwwave_acqus2xps(data_path, xps_temp_fn);

s.nii_fn = nii_temp_fn;
load(xps_temp_fn);
s.xps = xps;

% Something wrong with multiple echo acquisition so subsample shortest TE
mdm_s_subsample(s, xps.te == min(xps.te));

nii_sub_fn = fullfile(out_path, ['data_temp_sub' opt.nii_ext]);
xps_sub_fn = fullfile(out_path, 'data_temp_sub_xps.mat');

nii_fn = fullfile(out_path, ['data' opt.nii_ext]);
xps_fn = fullfile(out_path, 'data_xps.mat');

movefile(nii_sub_fn,nii_fn);
movefile(xps_sub_fn,xps_fn);

delete(nii_temp_fn);
delete(xps_temp_fn);

