function s_topup = mdm_s_topup(s_ap, s_pa, wp, o_fn, opt)
% function s_topup = mdm_s_topup(s_ap, s_pa, o_fn, wp, opt)
% 
% s_ap,   data acquired with phase coding anterior-posterior
% s_pa,   opposite phase coding polarity (posterior-anterior)
% wp,     working path
% o_fn,   output filename (with .nii.gz)
% opt,    option structure

if (nargin < 5), opt = []; end

opt = mdm_opt(opt);

% make output directories
msf_mkdir(wp); 
msf_mkdir(fileparts(o_fn));

% check which data to include in final apply topup
pa_fn = s_pa.nii_fn;
ap_fn = s_ap.nii_fn;

if (s_ap.xps.n == s_pa.xps.n)
    if (~all(s_ap.xps.b == s_pa.xps.b))
        error('ap and pa too different to merge by applytopup');
    end
    
    xps = s_ap.xps;
elseif (s_ap.xps.n > s_pa.xps.n)
    pa_fn = [];
    xps = s_ap.xps;
elseif (s_ap.xps.n < s_pa.xps.n)
    ap_fn = [];
    xps = s_pa.xps;
end

% Pull out the b0 / low b0 volumes
s_pa_b0 = mdm_s_subsample(s_pa, s_pa.xps.b <= 0.11e9, wp, opt);
s_ap_b0 = mdm_s_subsample(s_ap, s_ap.xps.b <= 0.11e9, wp, opt);

% Prepare powder average and merge (average all)
s_pa_b0.xps.a_ind = s_pa_b0.xps.b >= 0; 
s_ap_b0.xps.a_ind = s_ap_b0.xps.b >= 0;

s_pa_b0.xps = msf_rmfield(s_pa_b0.xps, 's_ind');
s_ap_b0.xps = msf_rmfield(s_ap_b0.xps, 's_ind');

% Powder average and merge 
s_pa_b0 = mdm_s_powder_average(s_pa_b0, wp, opt);
s_ap_b0 = mdm_s_powder_average(s_ap_b0, wp, opt);

s = mdm_s_merge({s_pa_b0, s_ap_b0}, wp, 'topup', opt);

% Write topup specification file 
% File is OK with both pa/ap and ap/pa order
% and define topup data path
topup_spec_fn = mdm_txt_write({'0 -1 0 1', '0 1 0 1'}, ...
    fullfile(wp, 'topup.txt'), opt);
topup_data_path = fullfile(wp, 'topup_data');

% run and apply topup
fsl_topup(s.nii_fn, topup_spec_fn, topup_data_path, opt);
fsl_applytopup(pa_fn, ap_fn, topup_data_path, topup_spec_fn, o_fn, opt);

% package it up and save the xps for good measure
s_topup.nii_fn = o_fn;
s_topup.xps = xps;
mdm_xps_save(s_topup.xps, mdm_xps_fn_from_nii_fn(s_topup.nii_fn));



