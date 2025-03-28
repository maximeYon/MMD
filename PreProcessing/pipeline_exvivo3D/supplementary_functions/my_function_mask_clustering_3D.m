function NewMask = my_function_mask_clustering_3D(Mask_path)
%% my manual edit mask
% Mask_path = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\ON-81-mbti-pilot-3d\16\pdata_mdd\nii_xps';

% Load data
opt = mdm_mrtrix_opt;
[data,h] = mdm_nii_read([Mask_path filesep 'data.nii.gz']);
Im_sum = sum(data,4);

% reshape and mean normalize data for clustering
data_size = size(data);
data = reshape(data,size(data,1)*size(data,2)*size(data,3),size(data,4));
data = data./mean(data);

%% perfom clustering
Ncluster = 12; % 12 for 3D

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
imagesc(squeeze(kmap(Nslice1,:,:)))
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
[~,idx] = max(numPixels);
for ind_obj = idx
    BW(CC.PixelIdxList{ind_obj}) = 1;
end
NewMask = zeros(size(NewMask));
NewMask(logical(BW)) =1;

%% display mask
figure(3)
for dim3 = 1:size(Im_sum,3)
    imagesc(squeeze(Im_sum(:,:,dim3)).*squeeze(NewMask(:,:,dim3)));
    axis image
    colormap gray
    drawnow
    title(num2str(dim3))
    pause(0.1)
end

%% Save mask
mdm_nii_write(uint8(NewMask), [Mask_path filesep 'data_mask.nii.gz'], h);
end




