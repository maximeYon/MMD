function my_function_split_data(data_path)

%% Open nifti data
data = niftiread(data_path);
header = niftiinfo(data_path);

b0s = data(:,:,:,1:2);
dataUp = data(:,:,:,[1, 3:end]);

%% create path b0s
path_only = split(data_path,filesep);
path_only = join(path_only(1:end-1,1),filesep,1); path_only = path_only{1};

path_b0s = [path_only filesep 'b0s.nii.gz'];

%% re-save NIFTI
my_save_NIFTI(b0s,header,path_b0s);
my_save_NIFTI(dataUp,header,data_path);

end


function my_save_NIFTI(data,header,nii_fn)
% make nifti header
h = mdm_nii_h_empty;
sdim = size(data);
h.pixdim(1+(1:length(sdim))) = sdim;
h.pixdim(2:4) = header.PixelDimensions(1:3);
h.xyzt_units = 'SI';
h.dim_info = [1 2 3];
h.sform_code = 0;
mdm_nii_write(data, nii_fn, h, 0);
end