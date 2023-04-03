%% Analyse Test_bvalue protocol
clearvars;close all;clc;
% Add function to path
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% Path for data
% filename = '20220825_mouse_brain_fin6\15'; 
filename = 'clay_1\13'; 
sorted = 'Yes'; % Yes

% format file names
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
filename = split(filename,'\'); filename = join(filename,filesep,1); filename = filename{1};
filename = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename  filesep 'pdata_mdd' filesep 'nii_xps' filesep 'data.nii.gz'];

% Load & squeeze nifti data
data = squeeze(niftiread(filename));

%% retrieve b-value and sort them
filename = split(filename,filesep); filename= join(filename(1:end-1,1),filesep,1); filename=filename{1};
load([filename filesep 'data_xps.mat']);

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
clearvars SortOrder2;

%% Compare ROI data
Nslice = 1;

figure(1)
% set(gcf,'color','w', 'Position', [6.60000000000000,726.600000000000,2040.80000000000,323.400000000000]);
LinData = data;
% LinData = reshape(LinData,sum(mask(:)),size(data,4));
% LinData = sum(LinData,1);
LinData = LinData./max(LinData);

im_ind = 1:size(LinData,1);

subplot(3,1,1)
plot(im_ind,xps.b(SortOrder)./10^9,'*');
subplot(3,1,2)
plot(im_ind,xps.te(SortOrder),'*');
subplot(3,1,3)
plot(im_ind,LinData(SortOrder),'*');



% % Display images
% Nslice = 27;
% x_Fov = 0:9.5/(size(data,1)-1):9.5;
% y_Fov = 0:16/(size(data,2)-1):16;
% ind_image = 1:size(SortOrder,1);
% ind_image = SortOrder(ind_image);
% figure(2)
% set(gcf,'color','w', 'Position', [0.1300    0.1100    0.7750    0.7417]);
% filename = 'brain_fin6_allim_expno27.gif';
% for ind = 1:numel(ind_image)
%     imagesc(x_Fov,y_Fov,data(:,:,Nslice,ind_image(ind))')
%     axis image
%     set(gca, 'YDir','normal')
%     colormap gray
%     title({['b-value = ' num2str(xps.b(ind_image(ind))./10^9)], ['TR = ' num2str(xps.tr(ind_image(ind)))], ['TE = ' num2str(xps.te(ind_image(ind))*1000)], ['image number ' num2str(ind) '/' num2str(numel(ind_image))] })
%     pause(0.2)
%     
%     %% Gif creation
%       frame = getframe(gcf);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if Nslice == 1
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%       end
% 
% end
