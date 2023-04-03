%% Top UP processing NIFTI dataset via FSL
clearvars;

home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

addpath('supplementary_functions');

% Define full paths to the nifti files to be analyzed
data_path = 'human_sample_fresh1\36';
mdd_path = 'pdata_mdd';

%% options
Display = 1;

data_path = split(data_path,filesep); data_path = join(data_path,filesep,1); data_path = data_path{1};
nii_fns{1} = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep mdd_path filesep 'nii_xps'];

%% get parameters from Paravision
data_path_pv = split(nii_fns{1},filesep);
data_path_pv = join(data_path_pv(1:end-2),filesep,1);
data_path_pv = [data_path_pv{1} filesep];
isTopUp = ReadPV360Param(data_path_pv, 'TopUpYesNo');
if ~strcmp(isTopUp,'yes')
    disp('Not Top Up data');
    return
end

%% Load xps & select images for b0 map calculation
load([nii_fns{1} filesep 'data_xps.mat']);
b0_indice = (xps.b==min(xps.b));%sum(b0_indice)
longTR_indice = (xps.tr>=2);%sum(longTR_indice)
shortTE_indice = (xps.te<=(min(xps.te)*3)); %sum(shortTE_indice)
selected_indices = b0_indice.*longTR_indice.*shortTE_indice;

%% Open nifti data
dataUp = niftiread([nii_fns{1} filesep 'dataUp.nii.gz']);
dataDown = niftiread([nii_fns{1} filesep 'dataDown.nii.gz']);
% size(dataUp)

%% reshape to Read phase slice rep
% Assume Read is the largest number
permuteYes = 0;
if size(dataUp,2)>=size(dataUp,1)
    dataUp = permute(dataUp,[2 1 3 4]);
    dataDown = permute(dataDown,[2 1 3 4]);
    permuteYes = 1;
end

%% Flip reversed blip data and correct phase shift
dataDown = flip(dataDown,2);

% figure(1)
% imshow3Dfull(dataUp(:,:,:,13));
% figure(2)
% imshow3Dfull(dataDown(:,:,:,13));

% correct phase offset here
Phase1Offset = ReadPV360Param(data_path_pv, 'PVM_SPackArrPhase1Offset');
FoV = ReadPV360Param(data_path_pv, 'PVM_Fov');
Matrix = ReadPV360Param(data_path_pv, 'PVM_Matrix');
if Phase1Offset ~=0
    pixel_shift = zeros(size(Matrix));
    pixel_shift(1,2) = (2*Phase1Offset)/FoV(1,2)*Matrix(1,2); % times 2 since it is applied in the wrong direction in PV processing
    for ind = 1:size(dataDown,4)
        dataDown(:,:,:,ind)=abs(fineshift(dataDown(:,:,:,ind),pixel_shift));
    end
end

%% repmat if single slice
SpatDim = ReadPV360Param(data_path_pv, 'PVM_SpatDimEnum');
NSlices = ReadPV360Param(data_path_pv, 'NSlices');
if ~strcmp(SpatDim,'<3d>')
    if NSlices ==1
        dataUp = repmat(dataUp,1,1,6,1);
        dataDown = repmat(dataDown,1,1,6,1);
    end
end

%% Get first images as b0 to compute the field map
% datab0= cat(4,dataUp(:,:,:,1),dataDown(:,:,:,1));
% datab0 = zeros(size(dataUp,1),size(dataUp,2),size(dataUp,3),sum(selected_indices)*2);
datab0 = zeros(size(dataUp,1),size(dataUp,2),size(dataUp,3),2);
datab0(:,:,:,1:2:end) = sum(dataUp(:,:,:,(selected_indices==1)),4);
datab0(:,:,:,2:2:end) = sum(dataDown(:,:,:,(selected_indices==1)),4);
% size(datab0)

%% display b0 images to check the offset corrections
if strcmp(SpatDim,'<3d>')
    scliceN =round(size(datab0,3)/2);
else
    scliceN =1;
end

if Display ==1
    figure(1)
    subplot(1,2,1)
    imagesc(squeeze(datab0(:,:,scliceN,1)))
    axis image
    title('blip up');
    subplot(1,2,2)
    imagesc(squeeze(datab0(:,:,scliceN,2)))
    axis image
    title('blip down');
    colormap gray
    drawnow
end

%% re-save NIFTI
my_path_NIFTI = split(nii_fns{1},filesep);
my_path_NIFTI = join(my_path_NIFTI(1:end),filesep,1);
my_path_NIFTI = [my_path_NIFTI{1} filesep];

my_save_NIFTI(datab0,data_path_pv,[my_path_NIFTI 'RefImg.nii.gz']);
my_save_NIFTI(dataDown,data_path_pv,[my_path_NIFTI 'dataDownFliped.nii.gz']);
my_save_NIFTI(dataUp,data_path_pv,[my_path_NIFTI 'dataUpFliped.nii.gz']);
% info = niftiinfo([my_path_NIFTI 'RefImg.nii.gz']);
% info.ImageSize %% Check that slice number is OK

%% perform FSL Top Up maps field maps calculation
disp('Performing Field map calculation')
tic;
Nimages_b0map = size(datab0,4)/2;
[~, result] = my_mdm_FSL_TopUp_FieldMap([my_path_NIFTI 'RefImg.nii.gz'], [my_path_NIFTI 'RefImgCorr.nii.gz'],Nimages_b0map);
if isempty(result)
    disp('Success')
else
    disp(result)
end
toc

%% Open Corrected Ref nifti images
data_corrb0 = niftiread([my_path_NIFTI 'RefImgCorr.nii.gz']);
field_map = niftiread([my_path_NIFTI 'fsl_fieldcoef.nii.gz']);

if ~strcmp(SpatDim,'<3d>')
    if NSlices ==1
        data_corrb0 = data_corrb0(:,:,3,:);
        field_map = field_map(:,:,3,:);
    end
end
% size(field_map)
% size(data_corrb0)

% display
if Display ==1
    figure(2)
    subplot(2,2,1)
    imagesc(squeeze(datab0(:,:,scliceN,1)))
    axis image
    title('blip up');
    subplot(2,2,2)
    imagesc(squeeze(datab0(:,:,scliceN,2)))
    axis image
    title('blip down');
    colormap gray
    subplot(2,2,3)
    imagesc(flip(squeeze(field_map(:,:,scliceN)),1))
    axis image
    title('Field map');
    subplot(2,2,4)
    imagesc(squeeze(data_corrb0(:,:,scliceN,1)))
    axis image
    title('Corrected image');
    colormap gray
    drawnow
end

%% perform FSL Top Up calculation of the entire image serie
disp('Performing Top Up correction on the image serie')
tic;
[~, result] = my_mdm_FSL_TopUp_Apply([my_path_NIFTI 'dataUpFliped.nii.gz'], [my_path_NIFTI 'dataDownFliped.nii.gz'],[my_path_NIFTI 'fsl'], [my_path_NIFTI 'ImgCorrTopUp.nii.gz']);
if isempty(result)
    disp('Success')
else
    disp(result)
end
toc

%% Open Corrected nifti data
data_corr = niftiread([my_path_NIFTI 'ImgCorrTopUp.nii.gz']);
if ~strcmp(SpatDim,'<3d>')
    if NSlices ==1
        data_corr = data_corr(:,:,3,:);
        dataUp = dataUp(:,:,3,:);
        dataDown = dataDown(:,:,3,:);
    end
end

%% display Original and corrected images
if Display ==1
    figure(3)
    for nrep = 1:min([size(data_corr,4) 10])
        subplot(2,2,1)
        imagesc(squeeze(dataUp(:,:,scliceN,nrep)))
        axis image
        title(['blip up n°' num2str(nrep)]);
        subplot(2,2,2)
        imagesc(squeeze(dataDown(:,:,scliceN,nrep)))
        axis image
        title(['blip down n°' num2str(nrep)]);
        colormap gray
        subplot(2,2,3)
        imagesc(flip(squeeze(field_map(:,:,scliceN)),1))
        axis image
        title('Field map');
        subplot(2,2,4)
        imagesc(squeeze(data_corr(:,:,scliceN,nrep)))
        axis image
        title(['Corrected image n°' num2str(nrep)]);
        colormap gray
        drawnow
        pause(0.1)
    end
end

%% save final data as data.nii.gz
if permuteYes ==1
    data_corr = permute(data_corr,[2 1 3 4]);
end
my_save_NIFTI(data_corr,data_path_pv,[my_path_NIFTI 'data.nii.gz']);

%% delete temporary files
delete([my_path_NIFTI 'fsl.nii.gz']);
delete([my_path_NIFTI 'imgcorrtopup.nii.gz']);
delete([my_path_NIFTI 'RefImg.nii.gz']);
delete([my_path_NIFTI 'refimgcorr.nii.gz']);
delete([my_path_NIFTI 'RefImg.topup_log']);
delete([my_path_NIFTI 'dataUpFliped.nii.gz']);
delete([my_path_NIFTI 'dataDownFliped.nii.gz']);



function my_save_NIFTI(data,data_path_pv,nii_fn)
PVM_SpatResol=ReadPV360Param(data_path_pv, 'PVM_SpatResol') ;
PVM_SliceThick=ReadPV360Param(data_path_pv, 'PVM_SliceThick') ;
% make nifti headear
h = mdm_nii_h_empty;
sdim = size(data);
h.pixdim(1+(1:length(sdim))) = sdim;
if numel(PVM_SpatResol) == 3 % For MSME with phase encode in both second and third dimensions
    h.pixdim(2:4) = [PVM_SpatResol];
end
if numel(PVM_SpatResol) == 2
    h.pixdim(2:4) = [PVM_SpatResol PVM_SliceThick];
end
h.xyzt_units = 'SI';
mdm_nii_write(data, nii_fn, h, 0);
end
