function my_function_Rician_Corr(data_path,doubleSampling,ncoilsSoS)

%% Open nifti data
data = niftiread(data_path);
header = niftiinfo(data_path);

%% Open nifti noise
path_only = split(data_path,filesep);
filename = path_only(end,1);filename = filename{1};
path_only = join(path_only(1:end-1,1),filesep,1); path_only = path_only{1};

if strcmp(filename,'dataDen.nii.gz')
    sigma = niftiread([path_only filesep 'sigma.nii.gz']);
end

%% perform Rician Bias correction
%  https://nyu-diffusionmri.github.io/DESIGNER-v2//docs/designer/background/#rician-bias-correction
% run.command('mrcalc working.mif 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if dwirc.mif')
data = data.^2;
if strcmp(doubleSampling,'yes')
    sigma = sigma.*(sqrt(ncoilsSoS*2)); % number of coils with sum_of_square addition * sum of square of the two double sampling images
else
    sigma = sigma.*(sqrt(ncoilsSoS)); % number of coils with sum_of_square addition
end
sigma = sigma.^2;
data2 = data-repmat(sigma,1,1,1,size(data,4)); 
data2 = abs(sqrt(data2));
data2 = data2.*isfinite(data2);

% data =sqrt(data);
% figure()
% subplot(1,2,1)
% hist(data(:),10000)
% subplot(1,2,2)
% hist(data2(:),10000)
% linkaxes

%% re-save NIFTI
my_save_NIFTI(data2,header,data_path);

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