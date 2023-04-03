function s = jm_MK_setup_xps_and_merge( meas_struct, output_path )

% meas_struct
%     .nii_img = cell array containing the path to each series in .nii
%     .diff_array = cell array containing double arrays with b-value, 
%                   b-vector, and lambda_axial for each img in respective
%                   series (read log-files with importdata)
% Index in diff_array needs to correspond to index in nii_img
    
% Total number of measurements to be analyzed
meas = size(meas_struct.nii_img, 2);

% Storage cell array
s = cell(meas, 1);

% Form struct with nifti path and xps information
for index = 1 : 1 : meas
    s{ index }.nii_fn   = meas_struct.nii_img{ index };
    s{ index }.xps      = jm_xps_from_log( meas_struct.diff_array{ index });
end

% Merge datasets
opt                 = mdm_opt();
opt.do_overwrite    = 1;
%wp                  = fullfile( pwd, 'work' );
wp                  = output_path;
nii_name            = 'MERGED_NII';
s                   = mdm_s_merge(s, wp, nii_name,opt);

end

