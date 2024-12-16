function mdm_bruker_PV360_2dseq2nii(data_path, nii_fn, rps, procno)
% function mdm_bruker_diff_wave_sym_msme_2dseq2nii(data_path, nii_fn, rps)
%
% Converting MSME images to nifti
% Image resolution in field n.pixdim in nifti header
%
% data_path: folder where the Bruker ser file is located
% nii_fn: nifti file name (including complete path and extension)
% rps: image recon parameters structure

if nargin == 1, rps = []; end

data_path_pv = [data_path filesep];
in_path = fullfile(data_path,'pdata',num2str(procno));

%% open reconstructed .mat file if NUS acquistion or spectroscopic data
PVM_SpatResol=ReadPV360Param(data_path_pv, 'PVM_SpatResol') ;
isNUS = ReadPV360Param(data_path_pv, 'CSYesNo');
% isNUS = 'yes'; %% force to take matlab file
if strcmp(isNUS,'yes')
    mat_filename = dir([in_path filesep 'im*.mat']);
    load([in_path filesep mat_filename(1,1).name]);
        if exist('imFT','var');imLLR = imFT;end
        if exist('imLLR_combined','var');imLLR = imLLR_combined;end
    imLLR = squeeze(imLLR);
    data = permute(imLLR,[2 1 3 4]);
    data = reshape(data,size(data,1),size(data,2),size(data,3),1,size(data,4));
else
    if numel(PVM_SpatResol) == 0 % localized spectroscopy 
    PV_version = ReadPV360Param(data_path_pv, 'ACQ_sw_version');
    if contains(PV_version,'pv 6.')==1
    in_path_fid = [data_path_pv 'rawdata.job0'];   
%     in_path_fid = [data_path_pv 'fid'];
    end
    if contains(PV_version,'pv-360')==1
    in_path_fid = [data_path_pv 'rawdata.job0'];    
    end
    fidID = fopen(in_path_fid,'r');
    np = ReadPV360Param(data_path_pv, 'PVM_SpecMatrix');
    nf = ReadPV360Param(data_path_pv, 'PVM_Nrepetitions'); 
    nAv = ReadPV360Param(data_path_pv, 'PVM_NAverages');
    fids = fread(fidID,'int32');
    fids = reshape(fids,2,np,nAv,nf);
    fids = squeeze(sum(fids,3));
    fids = squeeze(fids(1,:,:,:)-1i*fids(2,:,:,:));
%     fids=circshift(fids,[-68 0]);
    data = real(fftshift(fft(fids)));
%      plot(data)
    else % imaging
    imageObj = ImageDataObject(in_path);
    data = imageObj.data;
    end
end

% imshow3Dfull(real(data(:,:,:,1,1)))

% DwNDiffExp=ReadPVParam(data_path_pv, 'DwNDiffExp') ;
% MultipleEchos=ReadPVParam(data_path_pv, 'MultipleEchos') ;
% if (~isempty(DwNDiffExp) && strcmp(MultipleEchos,'yes'))
%     data = reshape(data,[size(data,1) size(data,2) size(data,3) size(data,5) size(data,6)]);
% end

% Size in various dimensions needs to be verified. DT 20201223
sz = size(data);
if numel(sz) == 6
    data = reshape(data,[sz(1), sz(2), sz(5), sz(6)]);
elseif numel(sz) == 2
    data = reshape(data,[sz(1), sz(2), 1, 1]);
elseif numel(sz) == 7
    data = reshape(data,[sz(1), sz(2), sz(5), sz(6)*sz(7)]);
else
    data = reshape(data,[sz(1), sz(2), sz(3), sz(4)*sz(5)]);
end

RECO_image_type=ReadPV360Param(data_path_pv, 'RECO_image_type', [data_path_pv 'pdata' filesep num2str(procno) filesep]);
if ~isreal(data); RECO_image_type = 'complex'; end
if strfind(RECO_image_type,'complex')
    data_real = real(data);
    data_imag = imag(data);
    data_real = flipdim(data_real,2);
    data_real = flipdim(data_real,1);
    data_imag = flipdim(data_imag,2);
    data_imag = flipdim(data_imag,1);
    data = complex(data_real,data_imag);
else
    data = flipdim(data,2);
    data = flipdim(data,1);
end


% if rps.denoising == 1
%     oldsize = size(data);
%     data = squeeze(data);
%     % create the mask
%     mask = zeros(size(data,1),size(data,2),size(data,3));
%     mask((data(:,:,1))> max(max(data(:,:,1)))/10) =1;
%     se = strel('disk',3);
%     mask = imclose(mask,se);
%     se = strel('disk',1);
%     mask = imdilate(mask,se);
%     mask = logical(mask);
%
%     %denoising
%     [dataden, ~, ~] = denoise(data, [3, 3, 2], mask);
%
%     dataFinale = reshape(dataden,oldsize);
%     data = dataFinale;
% end

PVM_SpatResol=ReadPV360Param(data_path_pv, 'PVM_SpatResol') ;
PVM_SliceThick=ReadPV360Param(data_path_pv, 'PVM_SliceThick') ;

%% Peak area extraction for spectroscopy
if numel(PVM_SpatResol) == 0 % localized spectroscopy 
data = my_bruker_spectro_extract(fids,data_path_pv);
end

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
if numel(PVM_SpatResol) == 0 % localized spectroscopy 
   PVM_VoxArrSize=ReadPV360Param(data_path_pv, 'PVM_VoxArrSize') ;
   h.pixdim(2:4) = PVM_VoxArrSize;
end
h.xyzt_units = 'SI';
h.dim_info = [1 2 3];
h.sform_code = 0;

% figure(12), clf, im2d = squeeze(data(:,:,2,1)); imagesc(im2d), return

% write nifti image and header

if strfind(RECO_image_type,'complex')
    data_real = real(data);
    data_imag = imag(data);
%     imshow3Dfull(data_imag(:,:,:,1,1))
    mdm_nii_write(data_real, [nii_fn(1:end-7) '_real' nii_fn(end-6:end)], h, 0);
    mdm_nii_write(data_imag, [nii_fn(1:end-7) '_imag' nii_fn(end-6:end)], h, 0);
else
    mdm_nii_write(data, nii_fn, h, 0);
end

if numel(PVM_SpatResol) == 0 % localized spectroscopy 
%% Save mask
NewMask = ones(size(data,1),size(data,2),size(data,3));
mdm_nii_write(uint8(NewMask), [nii_fn(1:end-17) filesep 'data_mask.nii.gz'], h);
end


