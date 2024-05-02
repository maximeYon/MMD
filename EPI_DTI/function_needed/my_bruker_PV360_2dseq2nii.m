function my_bruker_PV360_2dseq2nii(data_path, nii_fn, procno)
% Converting EPI images to nifti
% Image resolution in field n.pixdim in nifti header
%
% data_path: folder where the Bruker ser file is located
% nii_fn: nifti file name (including complete path and extension)
% rps: image recon parameters structure

if nargin == 1, rps = []; end

data_path_pv = [data_path filesep];
in_path = fullfile(data_path,'pdata',num2str(procno));

%% open reconstructed images
imageObj = ImageDataObject(in_path);
data = imageObj.data;

% Size in various dimensions needs to be verified.
sz = size(data);
if numel(sz) == 7
    data = permute(data,[1 2 3 4 5 7 6]);
     data = reshape(data,[sz(1), sz(2), sz(5), sz(6)*sz(7)]);   
elseif numel(sz) == 6
    data = reshape(data,[sz(1), sz(2), sz(5), sz(6)]);
elseif numel(sz) == 2
    data = reshape(data,[sz(1), sz(2), 1, 1]);
else
    data = reshape(data,[sz(1), sz(2), sz(3), sz(4)*sz(5)*sz(6)*sz(7)]);
end
data = flipdim(data,2);
data = flipdim(data,1);

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

% write nifti image and header
    mdm_nii_write(data, nii_fn, h, 0);
end