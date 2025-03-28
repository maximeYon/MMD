%% Launch Pre-processing on MDD data
clearvars; close all; clc;

addpath('supplementary_functions')
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% Define path list
path_list = {'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\tn28T_250203_2_5mmIprobe_CortOrg_1_4\5'}; % ;...
% path_list = {'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\27';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\28';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\29';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\30';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\31';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\32';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\33';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\34';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\35';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\36';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\37';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\38';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\39';...
%              'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain_CS\40'}; % ;...

force_top_up = 1;
%% Pre-processing loop
for ind = 1:size(path_list,1)
    %% set current data path
    data_path = path_list{ind};
    % data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];
    DoubleSampling = ReadPV360Param([data_path,filesep],'PVM_EpiCombine');
    ncoilsSoS = ReadPV360Param([data_path,filesep],'PVM_EncNReceivers');
    
    %% Step 2, create nifti and xps
    my_function_step2_pv2nii_xps(data_path)
    % Update data path
    data_path_pv = data_path;
    data_path = [data_path filesep 'pdata_mdd' filesep 'nii_xps'];
    pause(1)
    
    %% check if top up
    is_topup = ReadPV360Param([data_path_pv filesep], 'TopUpYesNo');
    is_firstOnly = ReadPV360Param([data_path_pv filesep], 'TopUpFirstOnly');
    SpatDim = ReadPV360Param([data_path_pv filesep], 'PVM_SpatDimEnum');
    
    if strcmp(is_topup,'yes')==1
        if strcmp(is_firstOnly,'yes')==1
            
            %% Split the data
            copyfile([data_path filesep 'data.nii.gz'], [data_path filesep 'dataRaw.nii.gz'], "f");
            my_function_split_data([data_path filesep 'dataRaw.nii.gz'])
            
            %%  denoise
            data_path_python = split(data_path,'\'); data_path_python = join(data_path_python,'/',1); data_path_python = data_path_python{1};
            if strcmp(SpatDim,'<3d>')
                my_function_python_denoising_3D(data_path_python);
            else
                my_function_python_denoising_2D(data_path_python);
            end
            
            %% Rician bias correction
            my_function_Rician_Corr([data_path filesep 'dataDen.nii.gz'],DoubleSampling,ncoilsSoS)
            
            %% deGibbs denoised data
            opt = mdm_mrtrix_opt;
            [status1, result1] = my_mdm_mrtrix_degibbs([data_path filesep 'dataDen.nii.gz'], [data_path filesep 'data.nii.gz'], opt);
            
            %% Motion correction
            %         my_function_MOCO_NIFTI([data_path filesep 'data.nii.gz'],'auto')
            
            %% TopUp
            my_function_topUp_singlePhase(data_path);
            pause(1)
            
            
        else
            %% split data into  blip up /down and denoise
            data_path_python = split(data_path,'\'); data_path_python = join(data_path_python,'/',1); data_path_python = data_path_python{1};
            my_function_python_denoisingtotpup_3D(data_path_python);
            pause(1)
            
            %% deGibbs denoised data
            disp('Start deGibbs')
            opt = mdm_mrtrix_opt;
            [status1, result1] = my_mdm_mrtrix_degibbs([data_path filesep 'denoised' filesep 'dataUp.nii.gz'], [data_path filesep 'dataUp.nii.gz'], opt);
            [status2, result2] = my_mdm_mrtrix_degibbs([data_path filesep 'denoised' filesep 'dataDown.nii.gz'], [data_path filesep 'dataDown.nii.gz'], opt);
            pause(1)
            
            %% TopUp
            my_function_topUp_NIFTI(data_path);
            pause(1)
            
            %% clean nii_xps folder
            delete([data_path filesep 'fsl_fieldcoef.nii.gz']);
            delete([data_path filesep 'dataUp.nii.gz']);
            delete([data_path filesep 'dataDown.nii.gz']);
            delete([data_path filesep 'fsl_movpar']);
            delete([data_path filesep 'denoised' filesep 'dataUp.nii.gz']);
            delete([data_path filesep 'denoised' filesep 'dataDown.nii.gz']);
            rmdir([data_path filesep 'denoised'],'s');
        end
    else
        copyfile([data_path filesep 'data.nii.gz'], [data_path filesep 'dataRaw.nii.gz'], "f");
        %%  denoise
        data_path_python = split(data_path,'\'); data_path_python = join(data_path_python,'/',1); data_path_python = data_path_python{1};
        my_function_python_denoising_3D(data_path_python);
        
        %% Rician bias correction
        my_function_Rician_Corr([data_path filesep 'dataDen.nii.gz'],DoubleSampling,ncoilsSoS)
        
        %% deGibbs denoised data
        opt = mdm_mrtrix_opt;
        [status1, result1] = my_mdm_mrtrix_degibbs([data_path filesep 'dataDen.nii.gz'], [data_path filesep 'data.nii.gz'], opt);
        
        %% Motion correction
        %         my_function_MOCO_NIFTI([data_path filesep 'data.nii.gz'],'auto')
        
        if force_top_up ==1

            %% TopUp
            my_function_topUp_singlePhase_forced(data_path);
            
        end
    end
    
    disp(['End of pre-processing ' num2str(ind)])
end

% %% Second loop for Mask creation
% for ind = 1:size(path_list,1)
%     %% set current data path
%     data_path = path_list{ind};
% %     data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];
%     % Update data path
%     data_path = [data_path filesep 'pdata_mdd' filesep 'nii_xps'];
%     %% Mask
%     my_function_mask_clustering_3D(data_path);
%     close all;
%     pause(1)
% end



