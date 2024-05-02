%% my manual edit mask
clearvars; close all; clc;
Mask_path = 'C:\Users\User\Documents\dataRennes\20240424_105148_pig_brain2\38\pdata_mdd\nii_xps\data.nii.gz';

addpath(genpath('Nifti')); % require NIFTI toolbox
% https://fr.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

% Load Nifti data and header
data = niftiread(Mask_path);
info = niftiinfo(Mask_path);
Im_sum = sum(data,4);

% reshape and mean normalize data for clustering
data_size = size(data);
data = reshape(data,size(data,1)*size(data,2)*size(data,3),size(data,4));
data = data./mean(data);

%% perfom clustering
Ncluster = 8; % 12 for 3D

opts = statset('Display','final','MaxIter',5);
GMModel = fitgmdist(data,Ncluster,'regularization',0.1,'Replicates',1,'Options',opts);
kmap = cluster(GMModel,data);

%% reshape kmap
kmap = reshape(kmap,data_size(1,1),data_size(1,2),data_size(1,3));

Nslice3 = round(size(kmap,3)/2);if Nslice3==0; Nslice3=1; end
Nslice2 = round(size(kmap,2)/2);if Nslice2==0; Nslice2=1; end
Nslice1 = round(size(kmap,1)/2);if Nslice1==0; Nslice1=1; end

figure(1)
set(gcf,'color','w','Position',[573.0000   75.6667  493.3333  782.3333]);
subplot(3,6,1:5)
imagesc(kmap(:,:,Nslice3))
axis image
subplot(3,6,7:11)
imagesc(squeeze(kmap(Nslice1,:,:)))
axis image
subplot(3,6,13:17)
imagesc(squeeze(kmap(:,Nslice2,:)))
colormap(hsv(Ncluster))
axis image
subplot(3,6,[12 18])
allColors = hsv(Ncluster);
hold on;
for ind = 1:Ncluster
scatter([],[],1, allColors(ind,:), 'filled', 'DisplayName', num2str(ind));
end
hold off
xticks([]);yticks([]);
axis ('off');
legend();

%% ask user to define noise area
prompt = {'Enter noise areas number separated by spaces'};
dlgtitle = 'Cluestering';
answer = inputdlg(prompt,dlgtitle);
indices = str2double(regexp(answer{1,1},'\d*','Match'));

%% Create Mask with kmap indice correspondig to the outside
NewMask = ones(size(kmap));
for indind= 1:numel(indices)
    NewMask(kmap==indices(indind)) =0;
end

%% close the mask
NewMask = imopen(NewMask,strel("sphere",2));
NewMask = imclose(NewMask,strel("sphere",5));

%% Select a single 3D object only
CC = bwconncomp(logical(NewMask),26);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
for ind_obj = idx
    BW(CC.PixelIdxList{ind_obj}) = 1;
end
NewMask = zeros(size(NewMask));
NewMask(logical(BW)) =1;

NewMask = imclose(NewMask,strel("disk",2));

%% display mask
figure(3)
for dim3 = 1:size(Im_sum,3)
    imagesc(squeeze(Im_sum(:,:,dim3)).*squeeze(NewMask(:,:,dim3)));
    axis image
    colormap gray
    drawnow
    title(num2str(dim3))
    pause(0.5)
end

%% Save mask
Save_path = split(Mask_path,filesep);
Save_path = join(Save_path(1:end-1,1),filesep,1); Save_path = Save_path{1};
nii = make_nii(NewMask);
save_nii(nii,[Save_path filesep 'my_mask.nii.gz']);
info2 = niftiinfo([Save_path filesep 'my_mask.nii.gz']);
info2.Transform.T = info.Transform.T;
niftiwrite(NewMask,[Save_path filesep 'my_mask'],info2,'Compressed',1)

%% Save mask single slice
Nslice = 30;
NewMaskSlice = zeros(size(NewMask));
NewMaskSlice(:,:,Nslice)= NewMask(:,:,Nslice);
niftiwrite(NewMaskSlice,[Save_path filesep 'my_maskSlice' num2str(Nslice)],info2,'Compressed',1)

%%
NewMask = zeros(size(NewMask));
% NewMask(28,71,4) = 1;
% NewMask(27,69:70,4) = 1;
% NewMask(36,28,4) = 1;
% NewMask(55:56,26:27,4) = 1;
% NewMask(57:61,82:85,4) = 1;
% NewMask(32,69:71,4) = 1;


