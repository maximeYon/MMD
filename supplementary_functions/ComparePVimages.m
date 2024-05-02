%% Compare two  immages
clearvars; close all; clc;

%% Define paths
pathData1 = 'C:\Users\Administrateur\Mon_Drive\Matlab\ProcessingPV360\data\20211220_tumor2\50\pdata\5';
pathData2 = 'C:\Users\Administrateur\Mon_Drive\Matlab\ProcessingPV360\data\20211220_tumor2\50\pdata\4';


% Import reconstructed images
imageObj1 = ImageDataObject(pathData1);
img1 = squeeze(imageObj1.data);
% img1= img1./mean(img1(:));

imageObj2 = ImageDataObject(pathData2);
img2 = squeeze(imageObj2.data);
% img2= img2./mean(img2(:));


%% Display
figure()
subplot(1,3,1)
imagesc(img1(:,:,1))
colorbar
subplot(1,3,2)
imagesc(img2(:,:,1))
colorbar
subplot(1,3,3)
imagesc(img1(:,:,1)-img2(:,:,1))
colorbar