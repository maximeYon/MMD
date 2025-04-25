%% Launch Pre-processing on MDD data
clearvars; close all; clc;

addpath('supplementary_functions')
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

path_list = {'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\pig_brain\10'}; % ;...
data_path = path_list{1};
data_path = [data_path filesep 'pdata_mdd' filesep 'nii_xps'];

% [status, result] =my_function_ANTs_Register([data_path filesep 'data_noTopUp.nii.gz'], [data_path filesep 'data_noTopUpRegisteredANTs.nii.gz'],1);
my_function_MOCO_NIFTI_rigid([data_path filesep 'data_noTopUp.nii.gz'],'No')

data_img_moving = squeeze(niftiread([data_path filesep 'data_noTopUp.nii.gz']));

ref = data_img_moving(:,:,:,1);
Nexp = size(data_img_moving,4);
img_corr = zeros(size(data_img_moving));
tformA= zeros(4,4,Nexp);

wb = waitbar(0,'Motion correction in progress');
for ind = 1:Nexp

[tform,reg] = imregmoment(data_img_moving(:,:,:,ind),ref,MedianThresholdBitmap=true); % )
img_corr(:,:,:,ind) =reg;
tformA(:,:,ind) = tform.A;
waitbar((ind/Nexp),wb,'Motion correction in progress');
end
close(wb)

Nslice =36;
figure(2)
for ind_img = 1:size(data_img_moving,4)
    subplot(1,2,1)
    imshowpair(data_img_moving(:,:,Nslice,ind_img)',squeeze(data_img_moving(:,:,Nslice,ind_img))',"Scaling","joint")
    title(['Without registration' num2str(ind_img)])
    subplot(1,2,2)
    imshowpair(data_img_moving(:,:,Nslice,ind_img)',squeeze(img_corr(:,:,Nslice,ind_img))',"Scaling","joint")
    title(['With registration' num2str(ind_img)])
    drawnow
    pause(1)
end