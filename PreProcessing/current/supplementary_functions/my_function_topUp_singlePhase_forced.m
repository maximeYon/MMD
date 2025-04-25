function my_function_topUp_singlePhase_forced(data_path)
%% Top UP processing NIFTI dataset via FSL
% options
Display = 1;

%% get parameters from Paravision
data_path_pv = split(data_path,filesep);
data_path_pv = join(data_path_pv(1:end-2,1),filesep,1); data_path_pv = data_path_pv{1};
data_path_pv = [data_path_pv filesep];
SpatDim = ReadPV360Param(data_path_pv, 'PVM_SpatDimEnum');
NSlices = ReadPV360Param(data_path_pv, 'NSlices');

%% Open nifti data
dataUp = niftiread([data_path filesep 'data.nii.gz']);
% size(dataUp)


%% reshape to Read phase slice rep
% Assume Read is the largest number
permuteYes = 0;
if size(dataUp,2)>=size(dataUp,1)
    dataUp = permute(dataUp,[2 1 3 4]);
    permuteYes = 1;
end

%% repmat if single slice
if ~strcmp(SpatDim,'<3d>')
    if NSlices < 10
        dataUp = reshape(dataUp,size(dataUp,1),size(dataUp,2),1,size(dataUp,3),size(dataUp,4));
        dataUp = repmat(dataUp,1,1,6,1,1);
        dataUp = cat(4,dataUp(:,:,:,1,:),dataUp(:,:,:,1,:),dataUp,dataUp(:,:,:,end,:),dataUp(:,:,:,end,:)); % double first and last slices
        dataUp = reshape(dataUp,size(dataUp,1),size(dataUp,2),size(dataUp,3)*size(dataUp,4),size(dataUp,5));
    end
end

if strcmp(SpatDim,'<3d>')
    scliceN =round(size(dataUp,3)/2);
else
    scliceN =1;
end

%% re-save NIFTI
my_path_NIFTI = [data_path filesep];
my_save_NIFTI_fakePixsize(dataUp,data_path_pv,[my_path_NIFTI 'dataFakePixsize.nii.gz']);
% info = niftiinfo([my_path_NIFTI 'RefImg.nii.gz']);
% info.ImageSize %% Check that slice number is OK


%% Open Corrected Ref nifti images
field_map = niftiread([my_path_NIFTI 'fsl_fieldcoef.nii.gz']);

if ~strcmp(SpatDim,'<3d>')
    if NSlices < 10
        Pos_slice = 15:6:size(data_corrb0,3)-12;
        field_map = field_map(:,:,Pos_slice,:);
    end
end

%% perform FSL Top Up calculation of the entire image serie
disp('Performing Top Up correction on the image serie')
tic;
[~, result] = my_function_FSL_TopUp_Apply_SinglePhase([my_path_NIFTI 'dataFakePixsize.nii.gz'],[my_path_NIFTI 'fsl'], [my_path_NIFTI 'dataCorrTopUp_jac.nii.gz']);

if isempty(result)
    disp('Success')
else
    disp(result)
end
toc

%% Open Corrected nifti data
data_corr = niftiread([my_path_NIFTI 'dataCorrTopUp_jac.nii.gz']);
if ~strcmp(SpatDim,'<3d>')
    if NSlices < 10
        Pos_slice = 15:6:size(data_corr,3)-12;
        data_corr = data_corr(:,:,Pos_slice,:);
        dataUp = dataUp(:,:,Pos_slice,:);
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
        colormap gray
        subplot(2,2,4)
        imagesc(squeeze(field_map(:,:,scliceN)))
        axis image
        title('Field map');
        subplot(2,2,3)
        imagesc(squeeze(data_corr(:,:,scliceN,nrep)))
        axis image
        title(['Corrected image n°' num2str(nrep)]);
        colormap gray
        drawnow
        pause(1)
    end
end

%% save final data as data.nii.gz
if permuteYes ==1
    data_corr = permute(data_corr,[2 1 3 4]);
end
my_save_NIFTI(data_corr,data_path_pv,[my_path_NIFTI 'data.nii.gz']);

%% delete temporary files
delete([my_path_NIFTI 'fsl.nii.gz']);
delete([my_path_NIFTI 'dataCorrTopUp_jac.nii.gz']);
delete([my_path_NIFTI 'dataFakePixsize.nii.gz']);

close all;

end

function my_save_NIFTI_fakePixsize(data,data_path_pv,nii_fn)
PVM_SpatResol=ReadPV360Param(data_path_pv, 'PVM_SpatResol') ;
PVM_SliceThick=ReadPV360Param(data_path_pv, 'PVM_SliceThick') ;
SliceGap=ReadPV360Param(data_path_pv, 'PVM_SPackArrSliceGap') ;
PVM_SliceThick = PVM_SliceThick+SliceGap(1,1);
% make nifti header
h = mdm_nii_h_empty;
sdim = size(data);
h.pixdim(1+(1:length(sdim))) = sdim;
if numel(PVM_SpatResol) == 3 % For MSME with phase encode in both second and third dimensions
    h.pixdim(2:4) = [PVM_SpatResol];
end
if numel(PVM_SpatResol) == 2
    h.pixdim(2:4) = [PVM_SpatResol PVM_SliceThick];
end

%% make resolution human like
h.pixdim(1,2:4) = (h.pixdim(1,2:4)./min(h.pixdim(1,2:4)))*2;

h.xyzt_units = 'SI';
h.dim_info = [1 2 3];
h.sform_code = 0;
mdm_nii_write(data, nii_fn, h, 0);
end

function my_save_NIFTI(data,data_path_pv,nii_fn)
PVM_SpatResol=ReadPV360Param(data_path_pv, 'PVM_SpatResol') ;
PVM_SliceThick=ReadPV360Param(data_path_pv, 'PVM_SliceThick') ;
SliceGap=ReadPV360Param(data_path_pv, 'PVM_SPackArrSliceGap') ;
PVM_SliceThick = PVM_SliceThick+SliceGap(1,1);
% make nifti header
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
h.dim_info = [1 2 3];
h.sform_code = 0;
mdm_nii_write(data, nii_fn, h, 0);
end