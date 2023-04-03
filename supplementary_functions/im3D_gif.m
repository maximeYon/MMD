%% Display 3D images
clearvars; close all; clc;
sorted = 'Yes'; % Yes

filename1 = '20230118_hippo3D\9';
filename = ['hippo3D_expno' filename1(1,end:end) '.gif'];
FoV = [9.5 14 9.5];

% format file names
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
filename1 = split(filename1,'\'); filename1 = join(filename1,filesep,1); filename1 = filename1{1};
filenameMask = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1 filesep 'pdata_mdd' filesep 'nii_xps' filesep  'data_mask.nii.gz'];
filename1 = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename1 filesep 'pdata_mdd' filesep 'nii_xps' filesep  'data.nii'];

% Load & squeeze nifti data
data1 = squeeze(niftiread(filename1));
mask = squeeze(niftiread(filenameMask));

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
data1 = data1(:,:,:,SortOrder);

% imshow3Dfull(data1(:,:,:,1));

%% Display substraction
ind_image = 1:size(SortOrder,1);
ind_image = SortOrder(ind_image);
x_im = 0:FoV(1,1)/(size(data1,1)-1):FoV(1,1);
y_im =flip(0:FoV(1,2)/(size(data1,2)-1):FoV(1,2));
z_im =0:FoV(1,3)/(size(data1,3)-1):FoV(1,3);
Slice_number = 35;
Slice_number2 = 67;
Slice_number3 = 47;
clims = [0 max(data1(:))/2.5];
figure(1)
set(gcf,'color','k', 'Position', [744.0000   79.4000  262.6000  970.6000]);
for im_ind = 1:size(data1,4)
    
    subplot(21,1,1:4)
    imagesc(y_im,x_im,squeeze(data1(:,:,Slice_number,im_ind)))
    axis image
    set(gca, 'Ydir', 'Normal')
    xticks([0 2 4 6 8 10 12 14 16])
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    title({['b-value = ' num2str(xps.b(ind_image(im_ind))./10^9)], ['TR = ' num2str(xps.tr(ind_image(im_ind)))], ['TE = ' num2str(xps.te(ind_image(im_ind))*1000)], ['image number ' num2str(im_ind) '/' num2str(size(data1,4))]},'Color','w')

    
    subplot(21,1,5:11)
    imagesc(z_im,x_im,flip(squeeze(data1(:,Slice_number2,:,im_ind))',1))
    axis image
    set(gca, 'Ydir', 'Normal')
    set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
    set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16],'YColor','white');
    colormap gray
    
    subplot(21,1,12:21)
    imagesc(z_im,y_im,flip(squeeze(data1(Slice_number3,:,:,im_ind)),1))
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
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
end

