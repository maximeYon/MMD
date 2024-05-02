%% Compare two  immages
clearvars; close all; clc;

%% Define paths
folder_path = 'C:\Users\Administrateur\Mon_Drive\Matlab\ProcessingPV360\data\20220204_tumor3\';
expno = [44 45 48 49 50 52 54];
procno = [1 1 1 1 1 1 1 1 1];
exp_descr = {'no VFA','VFA prepEsp 0.008','VFA prepEsp 0.01','VFA prepEsp 0.006','VFA Esp 0.002','VFA Esp 0.002 min 50°','VFA Esp 0.002 min 40°'};

%% Load data
Nexp = size(expno,2);
for ind = 1:Nexp
    pathData = [folder_path num2str(expno(ind)) filesep 'pdata' filesep num2str(procno(ind))];
    imageObj = ImageDataObject(pathData);
    img{ind} = permute(squeeze(imageObj.data),[2 1 3 4]);
end

%% Display
figure()
for ind = 1:Nexp
%     img_display= sum(img{ind},3);
    img_display= img{ind};
    img_display= img_display(:,:,3);
    subplot(2,Nexp,ind)
    imagesc(squeeze(img_display));
    axis image
    colormap gray
    title(exp_descr{ind})
    sign = img_display(22:38,22:38);
    noise = [img_display(1:8,1:8) img_display(57:64,1:8) img_display(1:8,57:64) img_display(57:64,57:64)];
    SNR(ind) = median(sign(:))/std(noise(:));
    if ind == Nexp
    subplot(2,Nexp,ind+1:2*ind)
    plot(SNR)
    end
end

