function out_path = mdm_bruker_dt_rare2d_recon(data_path, out_path, rps, opt)
% function out_path = mdm_bruker_dt_rare2d_recon(data_path, out_path, rps, opt)
%
% Recon images acquired with Bruker pulse sequence DT_axderare2d.
%
% Works in progress

if (nargin < 5), opt.present = 1; end

opt = mdm_opt(opt);

% Fit models
ind_start = 2; % Exclude first image

% Book-keeping

msf_mkdir(fullfile(out_path)); 
xps_fn = fullfile(out_path, 'xps.mat');
nii_fn = fullfile(out_path, ['data' opt.nii_ext]);

% convert bruker acquistion parameters to xps
mdm_bruker_acqus2mat(data_path);
mdm_bruker_dt_axderare2d_acqus2xps(data_path, xps_fn);
mdm_bruker_dt_rare2d_ser2nii(data_path, nii_fn, rps);

s.nii_fn = nii_fn;
load(xps_fn);
s.xps = xps;

mdm_s_subsample(s, (1:s.xps.n) >= ind_start);

