%% Downsampling
close all; clearvars; clc;
addpath('NIFTI_tools')

% Load data
data = load_nii('C:\Users\User\Desktop\review\MRI\N55134_25um_fa_RAS.nii');
img = data.img;
img = permute(img,[2 1 3]);

% Transform to Kspace
kspace = FFTXSpace2KSpace(FFTXSpace2KSpace(FFTXSpace2KSpace(img,1),2),3);
img_reco = abs(FFTKSpace2XSpace(FFTKSpace2XSpace(FFTKSpace2XSpace(kspace,1),2),3));

% % Display kspace
% figure(1)
% imagesc(abs(squeeze(kspace(Nslice,:,:))))

% Downsample
resolution = [30 40 50 100];
FoV = [18.8 12.8 12.8];
Res_orig = FoV./size(img);
for ind = 1:size(resolution,2)
    Npoints_down = round((FoV./(resolution(1,ind)/1000))/2)*2;
    diff_down = size(img)-Npoints_down;
    kspace_down = zeros(size(kspace));
    kspace_down(diff_down(1,1)/2:end-diff_down(1,1)/2,diff_down(1,2)/2:end-diff_down(1,2)/2,diff_down(1,3)/2:end-diff_down(1,3)/2)= kspace(diff_down(1,1)/2:end-diff_down(1,1)/2,diff_down(1,2)/2:end-diff_down(1,2)/2,diff_down(1,3)/2:end-diff_down(1,3)/2);
    data_reco{ind}.img = abs(FFTKSpace2XSpace(FFTKSpace2XSpace(FFTKSpace2XSpace(kspace_down,1),2),3));
end


%% Display
Nslice = 310;
figure(1)
subplot(2,size(resolution,2)+1,1)
imagesc(squeeze(img_reco(Nslice,:,:)));
title('25 \mu m')
axis image

subplot(2,size(resolution,2)+1,1+size(resolution,2)+1)
imagesc(squeeze(img(Nslice,:,:))-squeeze(img_reco(Nslice,:,:)));
axis image

for ind = 1:size(resolution,2)
    subplot(2,size(resolution,2)+1,ind+1)
    imagesc(squeeze(data_reco{ind}.img(Nslice,:,:)));
    title([num2str(resolution(1,ind)) ' \mu m'])
    axis image
    
    subplot(2,size(resolution,2)+1,ind+1+size(resolution,2)+1)
    imagesc(squeeze(img_reco(Nslice,:,:))- squeeze(data_reco{ind}.img(Nslice,:,:)));
    axis image 
    linkaxes
end  
colormap gray
    
    
    
    