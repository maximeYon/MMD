%% Display 3 D images
clearvars; close all; clc;
sorted = 'Yes'; % Yes

% filename1 = '20230118_hippo3D\10'; title1 = 'CS 28 %, 150 \mum iso';
% filename2 = '20230118_hippo3D\12'; title2 = 'CS 28 %, 150 \mum iso';
filename1 = '20230220_Rat_Brain_3D\17'; title1 = 'SE-EPI top-up, 25h';
filename2 = '20230220_Rat_Brain_3D\17'; title2 = 'SE-EPI no top up, 12.5 h';
filename = ['Rat_Brain_3D_expno' filename1(1,end-1:end) '_' filename2(1,end-1:end) '.gif'];
filename3Dscroll = ['3Dscroll_Rat_Brain_3D_expno' filename1(1,end-1:end) '_' filename2(1,end-1:end) '.gif'];

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
filename2 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename2 filesep 'pdata_mdd' filesep 'nii_xps' filesep  'orig_beforeTopUp' filesep 'dataBlipUpdn.nii'];

% Load & squeeze nifti data
data1 = squeeze(niftiread(filename1));
data2 = squeeze(niftiread(filename2));
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
    SortOrder_1 = SortOrder(SortOrder2);
else
    SortOrder_1 = 1:xps.n;
end
xps1 = xps;

clearvars SortOrder SortOrder2
filename2 = split(filename2,filesep); filename2= join(filename2(1:end-1,1),filesep,1); filename2=filename2{1};
load([filename2 filesep 'data_xps.mat']);
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
    SortOrder_2 = SortOrder(SortOrder2);
else
    SortOrder_2 = 1:xps.n;
end
xps2 = xps;
data1 = data1(:,:,:,SortOrder_1);
data2 = data2(:,:,:,SortOrder_2);


%% Display 3D
ind_image1 = 1:size(SortOrder_1,1);
ind_image1 = SortOrder_1(ind_image1);
ind_image2 = 1:size(SortOrder_2,1);
ind_image2 = SortOrder_2(ind_image2);
x_im = 0:FoV(1,1)/(size(data1,1)-1):FoV(1,1);
y_im =flip(0:FoV(1,2)/(size(data1,2)-1):FoV(1,2));
z_im =0:FoV(1,3)/(size(data1,3)-1):FoV(1,3);
Slice_number = 35;
Slice_number2 = 67;
Slice_number3 = 47;

x_im2 = 0:FoV(1,1)/(size(data2,1)-1):FoV(1,1);
y_im2 =flip(0:FoV(1,2)/(size(data2,2)-1):FoV(1,2));
z_im2 =0:FoV(1,3)/(size(data2,3)-1):FoV(1,3);
Slice_number_2 = 35;
Slice_number2_2 = 67;
Slice_number3_2 = 47;

% clims = [0 max(data1(:))/2.5];
figure(1)
set(gcf,'color','k', 'Position', [744,49.8,580.2,1020.8]);
for im_ind = 1:size(data1,4)
    
    subplot(21,2,[1 3 5 7])
    imagesc(y_im,x_im,squeeze(data1(:,:,Slice_number,im_ind)))
    axis image
    set(gca, 'Ydir', 'Normal')
    xticks([0 2 4 6 8 10 12 14 16])
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    title({title1,['b-value = ' num2str(xps1.b(ind_image1(im_ind))./10^9)], ['TR = ' num2str(xps1.tr(ind_image1(im_ind)))], ['TE = ' num2str(xps1.te(ind_image1(im_ind))*1000)], ['image number ' num2str(im_ind) '/' num2str(size(data1,4))] },'Color','w')

    subplot(21,2,[2 4 6 8])
    imagesc(y_im2,x_im2,squeeze(data2(:,:,Slice_number_2,im_ind)))
    axis image
    set(gca, 'Ydir', 'Normal')
    xticks([0 2 4 6 8 10 12 14 16])
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    title({title2,['b-value = ' num2str(xps2.b(ind_image2(im_ind))./10^9)], ['TR = ' num2str(xps2.tr(ind_image2(im_ind)))], ['TE = ' num2str(xps2.te(ind_image2(im_ind))*1000)], ['image number ' num2str(im_ind) '/' num2str(size(data2,4))] },'Color','w')
    
    subplot(21,2,[9 11 13 15 17 19 21])
    imagesc(z_im,x_im,flip(squeeze(data1(:,Slice_number2,:,im_ind))',1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    subplot(21,2,[10 12 14 16 18 22])
    imagesc(z_im2,x_im2,flip(squeeze(data2(:,Slice_number2_2,:,im_ind))',1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    subplot(21,2,23:2:41)
    imagesc(z_im,y_im,flip(squeeze(data1(Slice_number3,:,:,im_ind)),1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    subplot(21,2,24:2:42)
    imagesc(z_im2,y_im2,flip(squeeze(data2(Slice_number3_2,:,:,im_ind)),1))
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
end

%% Display 3D
figure(2)
clims1 = [0 max(max(max(data1(:,:,:,1))))/1.2];
clims2 = [0 max(max(max(data2(:,:,:,1))))/1.2];
set(gcf,'color','k', 'Position', [829.0000  477.8000  579.2000  407.2000]);
for slice_ind = 1:size(data2,3)
    subplot(1,2,1)
    imagesc(x_im,y_im,flip(squeeze(data1(:,:,round(slice_ind/size(data2,3)*size(data1,3)),1))'),clims1 )
    axis image
    set(gca, 'Ydir', 'Normal')
    xticks([0 2 4 6 8 10 12 14 16])
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    title({title1,['slice n°' num2str(round(slice_ind/size(data2,3)*size(data1,3)))]},'Color','w')
    
    subplot(1,2,2)
    imagesc(x_im,y_im,flip(squeeze(data2(:,:,slice_ind,1))'),clims2 )
    axis image
    set(gca, 'Ydir', 'Normal')
    xticks([0 2 4 6 8 10 12 14 16])
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    title({title2,['slice n°' num2str(slice_ind)]},'Color','w')
    
     frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if slice_ind == 1
        imwrite(imind,cm,filename3Dscroll,'gif','DelayTime',2, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename3Dscroll,'gif','DelayTime',2,'WriteMode','append');
    end
    
end
