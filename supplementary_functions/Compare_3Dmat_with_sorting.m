%% Analyse Test_bvalue protocol
clearvars;
% Add function to path
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% Path for data
filename = '20220722_mouse_brain_fin4\9'; 
sorted = 'Yes'; % Yes
mat_name1 = 'imLLR_0.01Lamb_30iter_56.25pct_W5_5';
mat_name2 = 'imLLR_0.01Lamb_30iter_56.25pct_W5_5first';

% format file names
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
filename = split(filename,'\'); filename = join(filename,filesep,1); filename = filename{1};
pathname = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename  filesep 'pdata' filesep '1' filesep];
filename = [home_path filesep 'ProcessingPV360' filesep 'data' filesep filename  filesep 'pdata_mdd' filesep 'nii_xps' filesep 'data.nii.gz'];

% Load & squeeze nifti data
load([pathname mat_name1 '.mat']);
im1 = imLLR; weights_factor1 = weights_factor;
clearvars imLLR weights_factor;
load([pathname mat_name2 '.mat']);
im2 = imLLR; weights_factor2 = weights_factor;
clearvars imLLR weights_factor;

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
Nslice = 32;

figure(1)
subplot(1,3,1)
imagesc(sum(im1(:,:,Nslice,:),4));
colormap gray
h = imfreehand;
mask = createMask(h);
mask_3D = zeros(size(im1));
mask_3D(:,:,Nslice,:) = repmat(mask,1,1,1,size(im1,4));

LinData1 = im1(logical(mask_3D));
LinData1 = reshape(LinData1,sum(mask(:)),size(im1,4));
LinData1 = sum(LinData1,1);
LinData1 = LinData1./max(LinData1);

LinData2 = im2(logical(mask_3D));
LinData2 = reshape(LinData2,sum(mask(:)),size(im2,4));
LinData2 = sum(LinData2,1);
LinData2 = LinData2./max(LinData2);

subplot(1,3,2:3)
hold on
plot(xps.b(SortOrder),LinData1(SortOrder),'*b');
plot(xps.b(SortOrder),LinData2(SortOrder),'*r');

% %% Display images
% ind_image = 1:4:306;
% ind_image = SortOrder(ind_image);
% figure(2)
% for ind = 1:numel(ind_image)
%     subplot(ceil(numel(ind_image)/10),10,ind)
%     imagesc(data(:,:,Nslice,ind_image(ind)))
%     colormap gray
%     title({['b-value = ' num2str(xps.b(ind_image(ind))./10^9)], ['TR = ' num2str(xps.tr(ind_image(ind)))], ['TE = ' num2str(xps.te(ind_image(ind))*1000)]})
% end
