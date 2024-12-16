function out_path = mdm_bruker_msme2nii_xps(data_path, out_path, rps, opt)
% function out_path = mdm_bruker_spen2nii_xps(data_path, out_path, rps, opt)
%
%
% Works in progress

if (nargin < 5), opt.present = 1; end

opt = mdm_opt(opt);

opt = mdm_opt();

nii_temp_fn = fullfile(out_path, ['data_temp' opt.nii_ext]);
xps_temp_fn = fullfile(out_path, 'data_temp_xps.mat');

mdm_bruker_diff_wave_sym_msme_2dseq2nii(data_path, nii_temp_fn, rps);
% mdm_bruker_diff_wave_sym_msme_acqus2xps(data_path, xps_temp_fn, rps);
% mdm_bruker_diff_wave_sym_msme_test_acqus2xps(data_path, xps_temp_fn, rps);
mdm_bruker_my_ll_R1_R2_DTD_MSME_PV360_acqus2xps(data_path, xps_temp_fn, rps);

s.nii_fn = nii_temp_fn;
load(xps_temp_fn);
s.xps = xps;

nii_fn = fullfile(out_path, ['data' opt.nii_ext]);
xps_fn = fullfile(out_path, 'data_xps.mat');

movefile(xps_temp_fn,xps_fn,'f');
movefile(nii_temp_fn,nii_fn,'f');

end


