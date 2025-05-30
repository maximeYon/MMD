%% Launch Pre-processing on MDD data
clearvars; close all; clc;

addpath('supplementary_functions')
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% Define path list
path_list = {'C:\Users\User\Documents\Data_kuopio_invivo\MMD_invivo_final\ON-93\ses-invivo\20230520_143622_ON_93_3d_1_1\16'}; % ;...
         
%% Pre-processing loop
for ind = 1:size(path_list,1)
    %% set current data path
    data_path = path_list{ind};
%     data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];
    
    %% Step 2, create nifti and xps
    my_function_step2_pv2nii_xps(data_path)
    % Update data path
    data_path = [data_path filesep 'pdata_mdd' filesep 'nii_xps'];
    pause(1)
    
    %% split data into  blip up /down and denoise
    data_path_python = split(data_path,'\'); data_path_python = join(data_path_python,'/',1); data_path_python = data_path_python{1};
    my_function_python_denoisingtotpup(data_path_python);
    pause(1)
    
    %% deGibbs denoised data
    opt = mdm_mrtrix_opt;
    [status1, result1] = my_mdm_mrtrix_degibbs([data_path filesep 'denoised' filesep 'dataUp.nii.gz'], [data_path filesep 'dataUp.nii.gz'], opt);
    [status2, result2] = my_mdm_mrtrix_degibbs([data_path filesep 'denoised' filesep 'dataDown.nii.gz'], [data_path filesep 'dataDown.nii.gz'], opt);
    pause(1)
    
    %% Drift correction
    my_function_MOCO_NIFTI_topup(data_path,'auto'); % secon argument is mode : if auto, no question MOCO is performed if shift > 0.8 pix
    pause(1)

    %% TopUp
    my_function_topUp_NIFTI(data_path);
    pause(1)
    
%     %% clean nii_xps folder
%     delete([data_path filesep 'dataMOCO.nii.gz']);
%     delete([data_path filesep 'fsl_fieldcoef.nii.gz']);
%     delete([data_path filesep 'dataUp.nii.gz']);
%     delete([data_path filesep 'dataDown.nii.gz']);  
%     delete([data_path filesep 'fsl_movpar']); 
%     delete([data_path filesep 'denoised' filesep 'dataUp.nii.gz']);
%     delete([data_path filesep 'denoised' filesep 'dataDown.nii.gz']); 
%     rmdir([data_path filesep 'denoised'],'s');  
%    
%     
    disp(['End of pre-processing ' num2str(ind)])
end

%% Second loop for Mask creation
for ind = 1:size(path_list,1)
    %% set current data path
    data_path = path_list{ind};
%     data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];
    % Update data path
    data_path = [data_path filesep 'pdata_mdd' filesep 'nii_xps'];
    %% Mask 
    my_function_mask_clustering(data_path); 
    close all;
    pause(1)
end



