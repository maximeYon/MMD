function mdm_bruker_diff_wave_sym_msme_2dseq2nii(data_path, nii_fn, rps)
% function mdm_bruker_diff_wave_sym_msme_2dseq2nii(data_path, nii_fn, rps)
%
% Converting MSME images to nifti
% Image resolution in field n.pixdim in nifti header
%
% data_path: folder where the Bruker ser file is located
% nii_fn: nifti file name (including complete path and extension)
% rps: image recon parameters structure

if nargin == 1, rps = []; end

data_path_pv = [data_path '/'];
in_path = fullfile(data_path,'pdata','1');

imageObj = ImageDataObject(in_path);
data = imageObj.data;

DwNDiffExp=ReadPVParam(data_path_pv, 'DwNDiffExp') ;
MultipleEchos=ReadPVParam(data_path_pv, 'MultipleEchos') ;
if (~isempty(DwNDiffExp) && strcmp(MultipleEchos,'yes'))
    data = reshape(data,[size(data,1) size(data,2) size(data,3) size(data,5) size(data,6)]);
end

% Size in various dimensions needs to be verified. DT 20201223
sz = size(data);
if numel(sz) == 6
    data = reshape(data,[sz(1), sz(2), sz(5), sz(6)]);
elseif numel(sz) == 2
    data = reshape(data,[sz(1), sz(2), 1, 1]);
else
    data = reshape(data,[sz(1), sz(2), sz(3), sz(4)*sz(5)]);
end
data = flipdim(data,2);
data = flipdim(data,1);

if rps.denoising == 1
    oldsize = size(data);
    data = squeeze(data);
    % create the mask
    mask = zeros(size(data,1),size(data,2),size(data,3));
    mask((data(:,:,1))> max(max(data(:,:,1)))/10) =1;
    se = strel('disk',3);
    mask = imclose(mask,se);
    se = strel('disk',1);
    mask = imdilate(mask,se);
    mask = logical(mask);

    %denoising
    [dataden, ~, ~] = denoise(data, [3, 3, 2], mask);

    dataFinale = reshape(dataden,oldsize);
    data = dataFinale;
end

PVM_SpatResol=ReadPVParam(data_path_pv, 'PVM_SpatResol') ;
PVM_SliceThick=ReadPVParam(data_path_pv, 'PVM_SliceThick') ; 

% make nifti headear
h = mdm_nii_h_empty;
sdim = size(data);
h.pixdim(1+(1:length(sdim))) = sdim;
if numel(PVM_SpatResol) == 3 % For MSME with phase encode in both second and third dimensions
    h.pixdim(2:4) = [PVM_SpatResol];
else
    h.pixdim(2:4) = [PVM_SpatResol PVM_SliceThick];
end
h.xyzt_units = 'SI';

% figure(12), clf, im2d = squeeze(data(:,:,2,1)); imagesc(im2d), return

% write nifti image and header
mdm_nii_write(data, nii_fn, h, 0);


