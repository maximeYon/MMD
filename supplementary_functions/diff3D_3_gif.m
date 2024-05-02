%% Display 3 D images
clearvars; close all; clc;
sorted = 'Yes'; % Yes

filename1 = '20230207_hippo3D\7'; title1 = '2 kHz VFA 4 cart';
filename2 = '20230207_hippo3D\14';  title2 = '10 kHz VFA 0 CS 36 echoes';
filename3 = '20230207_hippo3D\15';  title3 = '10 kHz VFA 0 CS 128 echoes';
filename = ['hippo3D_expno' filename1(1,end:end) '_' filename2(1,end:end) '_' filename3(1,end:end) '.gif'];
FoV = [9.5 14 9.5];

% format file names
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
filename1 = split(filename1,'\'); filename1 = join(filename1,filesep,1); filename1 = filename1{1};
% filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1  filesep 'pdata_mdd' filesep 'nii_xps' filesep 'orig' filesep 'data.nii.gz'];
% filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1 filesep 'pdata_mdd' filesep 'nii_xps'  filesep 'denoised' filesep  'data.nii.gz'];
filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1 filesep 'pdata_mdd' filesep 'nii_xps' filesep  'data.nii.gz'];
% filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1 filesep 'MRTrix_pdata_mdd' filesep 'nii_xps' filesep  'data.nii.gz'];

filename2 = split(filename2,'\'); filename2 = join(filename2,filesep,1); filename2 = filename2{1};
filename2 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename2 filesep 'pdata_mdd' filesep 'nii_xps' filesep  'data.nii.gz'];

filename3 = split(filename3,'\'); filename3 = join(filename3,filesep,1); filename3 = filename3{1};
filename3 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename3 filesep 'pdata_mdd' filesep 'nii_xps' filesep  'data.nii.gz'];

% Load & squeeze nifti data
data1 = squeeze(niftiread(filename1));
data2 = flip(squeeze(niftiread(filename2)),3);
data3 = flip(squeeze(niftiread(filename3)),3);
% data1 = data1./mean(data1(:));

%% retrieve b-value and sort them
filename1 = split(filename1,filesep); filename1= join(filename1(1:end-1,1),filesep,1); filename1=filename1{1};
load([filename1 filesep 'data_xps.mat']);
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
xps1 = xps;
filename2 = split(filename2,filesep); filename2= join(filename2(1:end-1,1),filesep,1); filename2=filename2{1};
load([filename2 filesep 'data_xps.mat']);
xps2 = xps;
filename3 = split(filename3,filesep); filename3= join(filename3(1:end-1,1),filesep,1); filename3=filename3{1};
load([filename3 filesep 'data_xps.mat']);
xps3 = xps;
data1 = data1(:,:,:,SortOrder);
data2 = data2(:,:,:,SortOrder);
data3 = data3(:,:,:,SortOrder);


%% Display 3D
ind_image = 1:size(SortOrder,1);
ind_image = SortOrder(ind_image);
x_im = 0:FoV(1,1)/(size(data1,1)-1):FoV(1,1);
y_im =flip(0:FoV(1,2)/(size(data1,2)-1):FoV(1,2));
z_im =0:FoV(1,3)/(size(data1,3)-1):FoV(1,3);
Slice_number = 25;
Slice_number2 = 65;
Slice_number3 = 30;
clims = [0 max(data1(:))/2.5];
figure(1)
set(gcf,'color','k', 'Position', [744,49.8,865,1020.8]);
for im_ind = 1:size(data1,4)
    
    subplot(21,3,1:3:3*4)
    imagesc(y_im,x_im,squeeze(data1(:,:,Slice_number,im_ind)))
    axis image
    set(gca, 'Ydir', 'Normal')
    xticks([0 2 4 6 8 10 12 14 16])
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    title({title1, ['b-value = ' num2str(xps1.b(ind_image(im_ind))./10^9)], ['TR = ' num2str(xps1.tr(ind_image(im_ind)))], ['TE = ' num2str(xps1.te(ind_image(im_ind))*1000)], ['image number ' num2str(im_ind) '/' num2str(size(data1,4))] },'Color','w')

    subplot(21,3,(1:3:3*4)+1)
    imagesc(y_im,x_im,squeeze(data2(:,:,Slice_number,im_ind)))
    axis image
    set(gca, 'Ydir', 'Normal')
    xticks([0 2 4 6 8 10 12 14 16])
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    title({title2,['b-value = ' num2str(xps2.b(ind_image(im_ind))./10^9)], ['TR = ' num2str(xps2.tr(ind_image(im_ind)))], ['TE = ' num2str(xps2.te(ind_image(im_ind))*1000)], ['image number ' num2str(im_ind) '/' num2str(size(data2,4))] },'Color','w')
    
    subplot(21,3,(1:3:3*4)+2)
    imagesc(y_im,x_im,squeeze(data3(:,:,Slice_number,im_ind)))
    axis image
    set(gca, 'Ydir', 'Normal')
    xticks([0 2 4 6 8 10 12 14 16])
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    title({title3,['b-value = ' num2str(xps3.b(ind_image(im_ind))./10^9)], ['TR = ' num2str(xps3.tr(ind_image(im_ind)))], ['TE = ' num2str(xps3.te(ind_image(im_ind))*1000)], ['image number ' num2str(im_ind) '/' num2str(size(data3,4))] },'Color','w')
    
    
    
    subplot(21,3,(13:3:13+3*7-3))
    imagesc(z_im,x_im,flip(squeeze(data1(:,Slice_number2,:,im_ind))',1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    subplot(21,3,(13:3:13+3*7-3)+1)
    imagesc(z_im,x_im,flip(squeeze(data2(:,Slice_number2,:,im_ind))',1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    subplot(21,3,(13:3:13+3*7-3)+2)
    imagesc(z_im,x_im,flip(squeeze(data3(:,Slice_number2,:,im_ind))',1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    
    
    
    subplot(21,3,(34:3:34+3*10-3))
    imagesc(z_im,y_im,flip(squeeze(data1(Slice_number3,:,:,im_ind)),1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    subplot(21,3,(34:3:34+3*10-3)+1)
    imagesc(z_im,y_im,flip(squeeze(data2(Slice_number3,:,:,im_ind)),1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    subplot(21,3,(34:3:34+3*10-3)+2)
    imagesc(z_im,y_im,flip(squeeze(data3(Slice_number3,:,:,im_ind)),1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    drawnow
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if im_ind == 1
        imwrite(imind,cm,filename,'gif','DelayTime',2, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',2,'WriteMode','append');
    end
   pause(2)
    
end

