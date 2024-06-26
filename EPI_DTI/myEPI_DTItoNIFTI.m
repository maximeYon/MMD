% Convert Bruker ParaVision data to nii and xps
clearvars;close all; clc;

% Add function to path
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));
addpath(genpath('function_needed'));

%% data path
data_path_file = 'ON-81-mbti-pilot-3d\17';
procno =1;

data_path = split(data_path_file,'\'); data_path = join(data_path,filesep,1); data_path = data_path{1};
data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];

if ~ isfile([data_path filesep 'dti.mat'])
    %% Convert data to NIFTI
    nii_fn = [data_path filesep 'dataUpDown.nii.gz'];
    my_bruker_PV360_2dseq2nii(data_path, nii_fn,procno);
    
    %% denoise data and split into up and down
    data_path_file = split(data_path_file,'\'); data_path_file = join(data_path_file,'/',1); data_path_file = data_path_file{1};
    my_launch_python_denoisingEPItotpup(data_path_file);
    
    %% deGibbs denoised data
    opt = mdm_mrtrix_opt;
    [status1, result1] = my_mdm_mrtrix_degibbs([data_path filesep 'denoised' filesep 'dataUp.nii.gz'], [data_path filesep 'dataUp.nii.gz'], opt);
    [status2, result2] = my_mdm_mrtrix_degibbs([data_path filesep 'denoised' filesep 'dataDown.nii.gz'], [data_path filesep 'dataDown.nii.gz'], opt);
    
    %% TopUp
    my_topUp_EPI_NIFTI(data_path);
    
    %% compute DTI
    [ADC,FA,reD,greeN,bluE,~,~] = DTI_EPI_PV_NIFTI(data_path);
    
    %% create mask
    mask = my_create_mask_clustering(data_path);
    
    %% Save DTI data
   save([data_path filesep 'dti.mat'],'ADC','FA','reD','greeN','bluE','data_path','mask');
else
    load([data_path filesep 'dti.mat'])
end

%% Display
n_slice = size(ADC,3);
RGB = zeros([3 size(ADC)]);
ADC = abs(ADC).*mask; FA = abs(FA).*mask;
RGB(1,:,:,:) = FA.*reD; RGB(2,:,:,:) = FA.*bluE; RGB(3,:,:,:) = FA.*greeN;
Fov=ReadPV360Param([data_path filesep], 'PVM_Fov') ;
y = 0:Fov(1,1)/(size(ADC,2)-1):Fov(1,1);
x = 0:Fov(1,2)/(size(ADC,1)-1):Fov(1,2);
clim_Diso = [0 0.001];
clim_FA = [0 1];
figure(100)
set(gcf,'color','k','Position',[663.4,49.8,1274.4,1020.8])
for ind_s = 1:n_slice
    subplot(3,size(ADC,3),ind_s)
    imagesc(x,y,flip(abs(squeeze(ADC(:,:,ind_s)))',2),clim_Diso)
    axis image
    set(gca, 'YDir', 'normal')
    set(gca,'xcolor', 'w'); set(gca,'ycolor', 'w');
    colormap gray
    
    subplot(3,size(FA,3),ind_s+n_slice)
    imagesc(x,y,flip(abs(squeeze(FA(:,:,ind_s)))',2),clim_FA)
    axis image
    set(gca, 'YDir', 'normal')
    set(gca,'xcolor', 'w'); set(gca,'ycolor', 'w');
    colormap gray
    
    subplot(3,size(FA,3),ind_s+n_slice*2)
    imagesc(x,y,flip(permute(abs(squeeze(RGB(:,:,:,ind_s))),[3 2 1]),2))
    axis image
    set(gca, 'YDir', 'normal')
%     set(gca, 'box','off');
    set(gca,'xcolor', 'w'); set(gca,'ycolor', 'w');
end

export_fig([data_path filesep 'result_dti'],'-bmp','-m3');
