function NewMask = my_function_mask_clustering(Mask_path)
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
Ncluster = 7; % 12 for 3D


opts = statset('Display','final','MaxIter',200);
GMModel = fitgmdist(data,Ncluster,'regularization',0.01,'Replicates',1,'Options',opts);
kmap = cluster(GMModel,data);

%% reshape kmap
kmap = reshape(kmap,data_size(1,1),data_size(1,2),data_size(1,3));

Nslice3 = round(size(kmap,3)/2);if Nslice3==0; Nslice3=1; end
Nslice1 = round(size(kmap,1)/2);if Nslice1==0; Nslice1=1; end
figure(20)
subplot(2,5,1:5)
imagesc(kmap(:,:,Nslice3))
axis image
subplot(2,5,6:9)
imagesc(squeeze(kmap(Nslice1,:,:)))
colormap(jet(Ncluster))
axis image
subplot(2,5,10)
allColors = jet(Ncluster);
hold on;
for ind = 1:Ncluster
scatter([],[],1, allColors(ind,:), 'filled', 'DisplayName', num2str(ind));
end
hold off
legend();


%% ask user to define noise area
prompt = {'Enter noise zone'};
dlgtitle = 'Cluestering';
answer = inputdlg(prompt,dlgtitle);
indices = str2double(regexp(answer{1,1},'\d*','Match'));


%% Create Mask with kmap indice correspondig to the outside
NewMask = ones(size(kmap));
% indices = [1 5 7];
for indind= 1:numel(indices)
    NewMask(kmap==indices(indind)) =0;
end

%% close the mask
NewMask = imopen(NewMask,strel("sphere",2));
NewMask = imclose(NewMask,strel("sphere",5));

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
mdm_nii_write(uint8(NewMask), [Mask_path filesep 'data_mask.nii.gz'], h);
end




