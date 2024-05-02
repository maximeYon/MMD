%% Compare two  immages
clearvars; close all; clc;

%% Define paths
folder_path = 'C:\Users\Administrateur\Mon_Drive\Matlab\ProcessingPV360\data\20220204_tumor3\';
expno = [26 27 29 28];
procno = [1 1 1 1];
exp_descr = {'VFA hard pulse','classical hard pulse','VFA BIR-4','classical BIR-4'};
slice = 16;
FoV = [20 9];


%% Load data
Nexp = size(expno,2);
for ind = 1:Nexp
    pathData = [folder_path num2str(expno(ind)) filesep 'pdata' filesep num2str(procno(ind))];
    imageObj = ImageDataObject(pathData);
    img{ind} = permute(squeeze(imageObj.data),[2 1 3 4]);
end

%% Display
y = 0:FoV(1,1)/(size(img{1},1)-1):FoV(1,1);
x = 0:FoV(1,2)/(size(img{1},2)-1):FoV(1,2);
figure()
% for rep = 1:166
    for ind = 1:Nexp
        img_display= sum(img{ind},4);
%         img_diplay= img{ind};
        img_display= squeeze(img_display(:,:,slice));
        subplot(2,2,ind)
        imagesc(x,y,squeeze(img_display));
        axis image
        colormap gray
        title(exp_descr{ind})
    end
%     drawnow
%     pause(1)
% end
