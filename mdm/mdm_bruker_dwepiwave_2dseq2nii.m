function res = mdm_bruker_dwepiwave_2dseq2nii(data_path, nii_fn, rps)
% function res = mdm_bruker_dwepiwave_2dseq2nii(data_path, nii_fn, rps)
%
% Converting EPI images to nifti
% Image resolution in field n.pixdim in nifti header
% Relies on Bruker's pvmatlab toolbox
%
% data_path: folder where the Bruker ser file is located
% nii_fn: nifti file name (including complete path and extension)
% rps: image recon parameters structure

if nargin < 3, rps.denoising = 0; end

if ~any(strcmp(mdm_bruker_readpvparam(data_path, 'PULPROG'),{lower('<rFOV_DWEpiWavev1_04.ppg>'); lower('<mcw_DWEpiWavev7.ppg>')})), res = 0; return, end
if ~strcmp(mdm_bruker_readpvparam(data_path, 'ACQ_pipe_status'),lower('Wrapup')), res = 0; return, end

in_path = fullfile(data_path,'pdata','1');

try
    imageObj = ImageDataObject(in_path); %Script from Bruker's pvmatlab toolbox
catch
    res = 0;
    warning('Conversion of image data failed. Check that the Bruker pvmatlab toolbox is in the Matlab serach path!')
    return
end

data = imageObj.data;

sz = size(data);
if numel(sz) == 5
    data = reshape(data,[sz(1), sz(2), sz(3), sz(5)]);
elseif numel(sz) == 6
    data = reshape(data,[sz(1), sz(2), sz(5), sz(6)]);
end
data = flipdim(data,2); data = flipdim(data,1);

if rps.denoising == 1
    oldsize = size(data);
    data = squeeze(data);
    % create the mask
    mask = zeros(size(data,1),size(data,2));
    mask((data(:,:,1))> max(max(data(:,:,1)))/10) =1;
    se = strel('disk',3);
    mask = imclose(mask,se);
    se = strel('disk',1);
    mask = imdilate(mask,se);
    mask = logical(mask);

    %denoising
    [dataden, ~, ~] = denoise(data, [2, 2], mask);

    dataFinale = reshape(dataden,oldsize);
    data = dataFinale;
end

PVM_SpatResol=mdm_bruker_readpvparam(data_path, 'PVM_SpatResol') ;
PVM_SliceThick=mdm_bruker_readpvparam(data_path, 'PVM_SliceThick') ;

% make nifti headear
h = mdm_nii_h_empty;
sdim = size(data);
h.pixdim(1+(1:length(sdim))) = sdim;
h.pixdim(2:4) = [PVM_SpatResol PVM_SliceThick];
h.xyzt_units = 'SI';
h.datatype = 16; %single precision

% write nifti image and header
if (~isempty(nii_fn))
    msf_mkdir(fileparts(nii_fn));
    mdm_nii_write(data, nii_fn, h, 0);  
end

res = 1;

