%% Display 3D images
clearvars; close all; clc;
sorted = 'Yes'; % Yes

filenames = {'EPItopup\53';'EPItopup\52'};
my_titles = {'2 segments';'Top Up'};
filenameGif = 'segVStopup.gif';
FoV = [15 9.5];

% format file names
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};

% Load & squeeze nifti data
Nimages = size(filenames,1);
for ind_im = 1:Nimages
%     %% Nifti
    filename_tmp = split(filenames{ind_im,1},'\'); filename_tmp = join(filename_tmp,filesep,1); filename_tmp = filename_tmp{1};
    filename_tmp = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename_tmp filesep 'pdata_mdd' filesep 'nii_xps' filesep  'data.nii.gz'];
    data{ind_im} = permute(squeeze(niftiread(filename_tmp)),[2 1 3 4]);
    %% Paravision
%     filename_tmp = split(filenames{ind_im,1},'\'); filename_tmp = join(filename_tmp,filesep,1); filename_tmp = filename_tmp{1};
%     filename_tmp = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename_tmp filesep 'pdata' filesep '1' filesep];
%     imageObj = ImageDataObject(filename_tmp);
%     data{ind_im} = permute(squeeze(imageObj.data),[2 1 3 4]);
end

%% retrieve b-value and sort them, possible only if NIFTI
filename_tmp = split(filename_tmp,filesep); filename_tmp= join(filename_tmp(1:end-1,1),filesep,1); filename_tmp=filename_tmp{1};
load([filename_tmp filesep 'data_xps.mat']);
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
for ind_im = 1:Nimages
data{ind_im} = data{ind_im}(:,:,SortOrder);
end


%% Display images
SliceN = 1;
for ind_im = 1:Nimages
    data_tmp = data{ind_im}(:,:,:,SliceN);
    clims(:,ind_im) = [0 median(data_tmp(:))].*40;
end
x_im = 0:FoV(1,1)/(size(data{1,1},1)-1):FoV(1,1);
y_im =flip(0:FoV(1,2)/(size(data{1,1},2)-1):FoV(1,2));
figure(1)
set(gcf,'color','k');
for ind_rep = 1; %:size(data{1,1},3)
    for ind_im = 1:Nimages
        subplot(1,Nimages,ind_im)
        imagesc(y_im,x_im,squeeze(data{ind_im}(:,:,ind_rep,SliceN)),clims(:,ind_im)')
        axis image
        set(gca, 'Ydir', 'Normal')
        xticks([0 2 4 6 8 10 12 14 16])
        yticks([0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30])
        set(gca,'XTickLabel',[0 2 4 6 8 10 12 14 16],'XColor','white');
        set(gca,'YTickLabel',[0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30],'YColor','white');
        colormap gray
        title({my_titles{ind_im},['image number ' num2str(ind_rep) '/' num2str(size(data{1,1},3))]},'Color','w')
%         title({my_titles{ind_im},['b-value = ' num2str(xps.b(ind_image(ind_rep))./10^9)], ['TR = ' num2str(xps.tr(ind_image(ind_rep)))], ['TE = ' num2str(xps.te(ind_image(ind_rep))*1000)], ['image number ' num2str(ind_rep) '/' num2str(size(data1,4))]},'Color','w')
    end
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ind_rep == 1
        imwrite(imind,cm,filenameGif,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filenameGif,'gif','WriteMode','append');
    end

end

