function [status, result] = my_mdm_Anima_TopUp_FieldMap(b0Up, b0Down)
% https://anima.readthedocs.io/en/latest/registration.html#epi-artifacts-correction
% 
% This function requires that Anima is installed on your computer. 

% animaDistortionCorrection -f AP_Image.nii.gz -b PA_Image.nii.gz -o Init_Correction.nii.gz -s 2

%% First Init calculation
path_data = [lower(b0Up(1)) strrep(b0Up(3:end-11),filesep,'/')];
reversed_dir = 0;

if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "animaDistortionCorrection -f ';
    cmd = [cmd '/mnt/' lower(b0Up(1)) strrep(b0Up(3:end),filesep,'/') '']; % input first b0Up
    cmd = [cmd ' -b /mnt/' lower(b0Down(1)) strrep(b0Down(3:end),filesep,'/') '']; % input then b0Down
    cmd = [cmd ' -o /mnt/' path_data 'Init_Correction.nii.gz -s 2 -d ' num2str(reversed_dir)];
else
cmd = 'animaDistortionCorrection -f ';
cmd = [cmd ' "' b0Up '"']; % input
cmd = [cmd ' -b "' b0Down '"']; % input
cmd = [cmd ' -o ' path_data 'Init_Correction.nii.gz -s 2 -d ' num2str(reversed_dir)];
end

if ispc == 1 % Suppose a WSL installation
    cmd = [cmd '"'];
    cmd = string(cmd);
end

[status, result] = msf_system(cmd);

% Alternative, copy to clipboard
% reduced_cmd = cmd{1}(13:end-1);
% clipboard('copy',reduced_cmd)

%% Second final calculation
% animaBMDistortionCorrection -f AP_Image.nii.gz -b PA_Image.nii.gz -i Init_Correction.nii.gz -o BM_Corrected_Image.nii.gz -O BM_Correction.nii.gz

if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "animaBMDistortionCorrection -f ';
    cmd = [cmd '/mnt/' lower(b0Up(1)) strrep(b0Up(3:end),filesep,'/') '']; % input first b0Up
    cmd = [cmd ' -b /mnt/' lower(b0Down(1)) strrep(b0Down(3:end),filesep,'/') '']; % input then b0Down
    cmd = [cmd ' -i /mnt/' path_data 'Init_Correction.nii.gz -o /mnt/' path_data 'BM_Corrected_Image.nii.gz -O /mnt/' path_data 'BM_Correction.nii.gz' ' -d ' num2str(reversed_dir)];
else
cmd = 'animaBMDistortionCorrection -f ';
cmd = [cmd ' "' b0Up '"']; % input
cmd = [cmd ' -b "' b0Down '"']; % input
cmd = [cmd ' -i ' path_data 'Init_Correction.nii.gz -o ' path_data 'BM_Corrected_Image.nii.gz -O ' path_data 'BM_Correction.nii.gz' ' -d ' num2str(reversed_dir)];
end

if ispc == 1 % Suppose a WSL installation
    cmd = [cmd '"'];
    cmd = string(cmd);
end

[status, result] = msf_system(cmd);

% Alternative, copy to clipboard
% reduced_cmd = cmd{1}(13:end-1);
% clipboard('copy',reduced_cmd)
end

