%% Compare two dataset directly
clearvars; close all; clc;
sorted = 'No'; % Yes

filename1 = '20221103_mouse_brain_fin7\20'; %classical
filename2 = '20221103_mouse_brain_fin7\20'; %NUS

% format file names
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
filename1 = split(filename1,'\'); filename1 = join(filename1,filesep,1); filename1 = filename1{1};
filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1  filesep 'pdata_mdd_full' filesep 'nii_xps' filesep 'data.nii.gz'];

filename2 = split(filename2,'\'); filename2 = join(filename2,filesep,1); filename2 = filename2{1};
filename2 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename2 filesep 'pdata_mdd_CS28' filesep 'nii_xps' filesep  'data.nii.gz'];

% Load & squeeze nifti data
data1 = squeeze(niftiread(filename1));
data2 = squeeze(niftiread(filename2));

%% retrieve b-value and sort them
filename2 = split(filename2,filesep); filename2= join(filename2(1:end-1,1),filesep,1); filename2=filename2{1};
load([filename2 filesep 'data_xps.mat']);
xpsEPI.b = xps.b;
xpsEPI.te = xps.te;
xpsEPI.tr = xps.tr;

filename1 = split(filename1,filesep); filename1= join(filename1(1:end-1,1),filesep,1); filename1=filename1{1};
load([filename1 filesep 'data_xps.mat']);
xpsRARE.b = xps.b;
xpsRARE.te = xps.te;
xpsRARE.tr = xps.tr;
if strcmp(sorted,'Yes')
    [~,SortOrder] = sort(xps.b);
    Bval_tmp = xps.b./1e9;
    Bval_tmp = Bval_tmp(SortOrder);
    Bdelta_tmp = xps.b_delta(SortOrder);
    BvalChange = find(diff(Bval_tmp)>0.01);
    for ind = 1:size(BvalChange,1)
        if ind ==1
           [~,SortOrder2(1:BvalChange(ind),:)] = sort(Bdelta_tmp(1:BvalChange(ind),:));
        else
            S_SortOrder2 = size(SortOrder2,1);
            [~,SortOrder2(BvalChange(ind-1)+1:BvalChange(ind),:)] = sort(Bdelta_tmp(BvalChange(ind-1)+1:BvalChange(ind),:));
            SortOrder2(BvalChange(ind-1)+1:BvalChange(ind),:) = SortOrder2(BvalChange(ind-1)+1:BvalChange(ind),:)+S_SortOrder2;
        end
        if ind ==size(BvalChange,1)
            S_SortOrder2 = size(SortOrder2,1);
           [~,SortOrder2(BvalChange(ind)+1:size(SortOrder,1),:)] = sort(Bdelta_tmp(BvalChange(ind)+1:end,:));
            SortOrder2(BvalChange(ind)+1:size(SortOrder,1),:) = SortOrder2(BvalChange(ind)+1:size(SortOrder,1),:)+S_SortOrder2;
        end
    end
    SortOrder = SortOrder(SortOrder2);
else
    SortOrder = 1:xps.n;
end

%% Normalization
% data1 = data1./mean(data1(:));
% data2 = data2./mean(data2(:));

imgSum1 = sum(data1,4);
imgSum2 = sum(data2,4);

imgTot = imgSum1+imgSum2;

%% Compare ROI data
figure(1)
Nslice = 25;
subplot(1,3,1)
imagesc(imgTot(:,:,Nslice))
colormap gray
h = imfreehand;
mask_slice = createMask(h);
mask = zeros(size(data1,1),size(data1,2),size(data1,3));
mask(:,:,Nslice) = mask_slice;
mask = logical(repmat(mask,1,1,1,size(data1,4)));

data1 = data1./max(data1(mask));
data2 = data2./max(data2(mask));

LinData1 = data1(mask);
LinData1 = reshape(LinData1,sum(mask(:))/size(data1,4),size(data1,4));
LinData1 = sum(LinData1,1);

LinData2 = data2(mask);
LinData2 = reshape(LinData2,sum(mask(:))/size(data2,4),size(data2,4));
LinData2 = sum(LinData2,1);

subplot(1,3,[2 3])
hold on
plot(SortOrder,LinData1(SortOrder),'*b')
plot(SortOrder,LinData2(SortOrder),'*r')
legend('classical','NUS');


%% Display substraction
figure(2)
if numel(size(data1))==4; dim3D =1; else dim3D =0; end
Slice_number = 32;
Slice_number2 = 32;
for im_ind = 1:size(data1,3)
figure(1)
subplot(1,3,1)
if dim3D ==1
imagesc(squeeze(data1(:,:,Slice_number,im_ind)))    
else
imagesc(squeeze(data1(:,:,im_ind)))
end
axis image
colormap gray
colorbar
title('orig')
subplot(1,3,2)
if dim3D ==1
imagesc(squeeze(data2(:,:,Slice_number2,im_ind)))    
else
imagesc(squeeze(data2(:,:,im_ind)))
end
axis image
colormap gray
colorbar
title('NUS 28%')
subplot(1,3,3)
if dim3D ==1
imagesc(squeeze(data1(:,:,Slice_number,im_ind))-squeeze(data2(:,:,Slice_number,im_ind))) 
else
imagesc(squeeze(data1(:,:,im_ind))-squeeze(data2(:,:,im_ind)))
end
axis image
colormap gray
colorbar
title('diff')
drawnow
pause(1)
end


