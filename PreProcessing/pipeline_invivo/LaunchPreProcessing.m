%% Launch Pre-processing on MDD data
clearvars; close all; clc;

addpath('supplementary_functions')
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% Define path list
path_list = {'invivo_methodo_Rician_corr_sqrt2_correct\6';}; %...
%              'invivo_Methodo_ON_83_28d_1_1\25';...
%              'invivo_Methodo_ON_83_28d_1_1\26';...
%              'invivo_Methodo_ON_83_28d_1_1\27'}; % ;...
         
%% Pre-processing loop
for ind = 1:size(path_list,1)
    %% set current data path
    data_path = path_list{ind};
    data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];
    % get parameters
    DoubleSampling = ReadPV360Param([data_path,filesep],'PVM_EpiCombine');
    ncoilsSoS = ReadPV360Param([data_path,filesep],'PVM_EncNReceivers');
    
    %% Step 2, create nifti and xps
    my_function_step2_pv2nii_xps(data_path)
    % Update data path
    data_path = [data_path filesep 'pdata_mdd' filesep 'nii_xps'];
    pause(1)
    
    %% split data into  blip up /down and denoise
    data_path_python = split(data_path,'\'); data_path_python = join(data_path_python,'/',1); data_path_python = data_path_python{1};
     my_function_python_denoisingtotpup(data_path_python);
% %         my_function_python_NOdenoisingtotpup(data_path_python);
    pause(1)
    
     %% Rician bias correction
    my_function_Rician_Corr([data_path filesep 'denoised' filesep 'dataUp.nii.gz'],DoubleSampling,ncoilsSoS)
    my_function_Rician_Corr([data_path filesep 'denoised' filesep 'dataDown.nii.gz'],DoubleSampling,ncoilsSoS)
    
    % %% deGibbs denoised data
    opt = mdm_mrtrix_opt;
    [status1, result1] = my_mdm_mrtrix_degibbs([data_path filesep 'denoised' filesep 'dataUp.nii.gz'], [data_path filesep 'dataUp.nii.gz'], opt);
    [status2, result2] = my_mdm_mrtrix_degibbs([data_path filesep 'denoised' filesep 'dataDown.nii.gz'], [data_path filesep 'dataDown.nii.gz'], opt);
    copyfile([data_path filesep 'denoised' filesep 'sigmaUp.nii.gz'],[data_path filesep 'sigmaUp.nii.gz'],'f')
    copyfile([data_path filesep 'denoised' filesep 'sigmaDown.nii.gz'],[data_path filesep 'sigmaDown.nii.gz'],'f')
    % pause(1)
    % 
    % %% Drift correction
    my_function_MOCO_NIFTI_topup(data_path,'auto'); % secon argument is mode : if auto, no question MOCO is performed if shift > 0.8 pix
    pause(1)

    %% TopUp
    my_function_topUp_NIFTI(data_path);
    pause(1)
    
    %% clean nii_xps folder
    delete([data_path filesep 'dataMOCO.nii.gz']);
    delete([data_path filesep 'fsl_fieldcoef.nii.gz']);
    delete([data_path filesep 'dataUp.nii.gz']);
    delete([data_path filesep 'dataDown.nii.gz']);  
    delete([data_path filesep 'fsl_movpar']); 
    delete([data_path filesep 'denoised' filesep 'dataUp.nii.gz']);
    delete([data_path filesep 'denoised' filesep 'dataDown.nii.gz']); 
    rmdir([data_path filesep 'denoised'],'s');  
   
    
    disp(['End of pre-processing ' num2str(ind)])
end

% %% Second loop for Mask creation
% for ind = 1:size(path_list,1)
%     %% set current data path
%     data_path = path_list{ind};
%     data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];
%     % Update data path
%     data_path = [data_path filesep 'pdata_mdd' filesep 'nii_xps'];
%     %% Mask 
%     my_function_mask_clustering(data_path); 
%     close all;
%     pause(1)
% end



