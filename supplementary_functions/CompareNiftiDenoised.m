%% Compare two dataset directly
clearvars; close all; clc;
sorted = 'Yes'; % Yes

filename1 = '20230118_hippo3D\19';
filename2 = '20230118_hippo3D\19';

% format file names
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
filename1 = split(filename1,'\'); filename1 = join(filename1,filesep,1); filename1 = filename1{1};
% filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1  filesep 'pdata_mdd' filesep 'nii_xps' filesep 'orig' filesep 'data.nii.gz'];
% filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1 filesep 'pdata_mdd' filesep 'nii_xps'  filesep 'denoised' filesep  'data.nii'];
filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1 filesep 'pdata_mdd' filesep 'nii_xps' filesep 'orig' filesep];
% filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1 filesep 'MRTrix_pdata_mdd' filesep 'nii_xps' filesep  'data.nii.gz'];

filename2 = split(filename2,'\'); filename2 = join(filename2,filesep,1); filename2 = filename2{1};
% filename2 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename2 filesep 'pdata_mdd' filesep 'nii_xps'  filesep 'denoised' filesep  'data.nii'];
filename2 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename2 filesep 'pdata_mdd' filesep 'nii_xps' filesep 'denoised' filesep  'data.nii'];

% Load & squeeze nifti data
data1_real = squeeze(niftiread([filename1 'data_real.nii.gz']));
data1_imag = squeeze(niftiread([filename1 'data_imag.nii.gz']));
data1 = abs(complex(data1_real,data1_imag));

data2 = squeeze(niftiread(filename2));

% data1 = data1./mean(data1(:));
% data2 = data2./mean(data2(:));

%% retrieve b-value and sort them
filename2 = split(filename2,filesep); filename2= join(filename2(1:end-1,1),filesep,1); filename2=filename2{1};
load([filename2 filesep 'data_xps.mat']);
xps1.b = xps.b;
xps1.te = xps.te;
xps1.tr = xps.tr;

filename1 = split(filename1,filesep); filename1= join(filename1(1:end-1,1),filesep,1); filename1=filename1{1};
load([filename1 filesep 'data_xps.mat']);
xps2.b = xps.b;
xps2.te = xps.te;
xps2.tr = xps.tr;
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


%% Display substraction
if numel(size(data1))==4; dim3D =1; else dim3D =0; end
Slice_number = 35;
Slice_number2 = 35;
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
title('denoised')
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

