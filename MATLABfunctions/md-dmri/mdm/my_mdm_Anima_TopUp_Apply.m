function [status, result] = my_mdm_Anima_TopUp_Apply(BlipUpImg, BlipDownImg)
% https://anima.readthedocs.io/en/latest/registration.html#epi-artifacts-correction
% 
% This function requires that Anima is installed on your computer. 

% animaApplyDistortionCorrection -f DWI_AP.nii.gz -b DWI_PA.nii.gz -t BM_Correction.nii.gz -o DWI_Corrected.nii.gz

path_data = [lower(BlipUpImg(1)) strrep(BlipUpImg(3:end-19),filesep,'/')];

if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "animaApplyDistortionCorrection -f ';
    cmd = [cmd '/mnt/' lower(BlipUpImg(1)) strrep(BlipUpImg(3:end),filesep,'/') '']; % input first b0Up
    cmd = [cmd ' -b /mnt/' lower(BlipDownImg(1)) strrep(BlipDownImg(3:end),filesep,'/') '']; % input then b0Down
    cmd = [cmd ' -t /mnt/' path_data 'BM_Correction.nii.gz -o /mnt/' path_data 'DWI_Corrected.nii.gz'];
else
cmd = 'animaApplyDistortionCorrection -f ';
cmd = [cmd ' "' b0Up '"']; % input
cmd = [cmd ' -b "' b0Down '"']; % input
cmd = [cmd ' -t ' path_data 'BM_Correction.nii.gz -o ' path_data 'DWI_Corrected.nii.gz'];
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

