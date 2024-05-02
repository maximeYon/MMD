function [ADC,FA,reD,greeN,bluE,VectorF,DifT] = DTI_EPI_PV_NIFTI_MRtrix(data_path)

%% Get parameters
data_path_pv = [data_path filesep];
SlicePerPack=ReadPV360Param(data_path_pv, 'PVM_SPackArrNSlices') ;
DwAoImages=ReadPV360Param(data_path_pv, 'PVM_DwAoImages') ;
DwNDiffDir=ReadPV360Param(data_path_pv, 'PVM_DwNDiffDir') ;
NDiffDir = DwNDiffDir+DwAoImages;
DwEffBval=ReadPV360Param(data_path_pv, 'PVM_DwEffBval') ;
DwGradVec=ReadPV360Param(data_path_pv, 'PVM_DwDir') ;
SpatDimEnum=ReadPV360Param(data_path_pv, 'PVM_SpatDimEnum') ;
Fov=ReadPV360Param(data_path_pv, 'PVM_Fov') ;
DwBMat=ReadPV360Param(data_path_pv, 'PVM_DwBMat') ;
%% Import reconstructed images
% data = niftiread([data_path filesep 'data.nii.gz']);

%% Prepare diff data
GradientOrientations = [zeros(DwAoImages,3);DwGradVec];
GO(:,1) = GradientOrientations(:,1);
GO(:,2) = GradientOrientations(:,2);
GO(:,3) = GradientOrientations(:,3);
GradientOrientations = GO;

%% Write b-value file and gradient orientation file
my_format = [];
for ind = 1:size(DwEffBval,2)-1
my_format = [my_format '%f '];
end
my_format = [my_format '%f\n'];
fileID = fopen([data_path filesep 'bvals.txt'],'w');
fprintf(fileID,my_format,DwEffBval)
fclose(fileID);

my_format = [];
for ind = 1:size(GradientOrientations,1)-1
my_format = [my_format '%f '];
end
my_format = [my_format '%f\n'];
fileID = fopen([data_path filesep 'bvecs.txt'],'w');
fprintf(fileID,my_format,GradientOrientations(:,2))
fprintf(fileID,my_format,GradientOrientations(:,1).*-1)
fprintf(fileID,my_format,GradientOrientations(:,3))
fclose(fileID);

%% Do DTI calculation
if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "dwi2tensor';
    cmd = [cmd ' /mnt/' lower(data_path(1)) strrep(data_path(3:end),'\','/') '/data.nii.gz']; % input
    cmd = [cmd ' /mnt/' lower(data_path(1)) strrep(data_path(3:end),'\','/') '/DTIout.nii.gz -force']; % output
    cmd = [cmd ' -fslgrad /mnt/' lower(data_path(1)) strrep(data_path(3:end),'\','/') '/bvecs.txt'];
    cmd = [cmd ' /mnt/' lower(data_path(1)) strrep(data_path(3:end),'\','/') '/bvals.txt'];
    cmd = [cmd '"'];
    cmd = string(cmd);
else
cmd = 'dwi2tensor';
cmd = [cmd ' "' data_path filesep 'data.nii.gz' '"']; % input
cmd = [cmd ' "' data_path filesep 'DTIout.nii.gz' '" -force']; % output
cmd = [cmd ' -fslgrad "' data_path filesep 'bvecs.txt' '"']; % output
cmd = [cmd ' "' data_path filesep 'bvals.txt' '"']; % output

end

[status, result] = msf_system(cmd);

%% Create metrics
if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "tensor2metric';
    cmd = [cmd ' -adc /mnt/' lower(data_path(1)) strrep(data_path(3:end),'\','/') '/adc.nii.gz']; % output
    cmd = [cmd ' -fa /mnt/' lower(data_path(1)) strrep(data_path(3:end),'\','/') '/fa.nii.gz']; % output
    cmd = [cmd ' -vector /mnt/' lower(data_path(1)) strrep(data_path(3:end),'\','/') '/vect.nii.gz -modulate none']; % output
    cmd = [cmd ' /mnt/' lower(data_path(1)) strrep(data_path(3:end),'\','/') '/DTIout.nii.gz -force']; % input
    cmd = [cmd '"'];
    cmd = string(cmd);
else
%cmd = 'tensor2metric -adc image -fa image -vector image';
cmd = 'tensor2metric';
cmd = [cmd ' -adc "' data_path filesep 'adc.nii.gz' '"']; % output
cmd = [cmd ' -fa "' data_path filesep 'fa.nii.gz' '"']; % output
cmd = [cmd ' -ad "' data_path filesep 'ad.nii.gz' '"']; % output
cmd = [cmd ' -rd "' data_path filesep 'rd.nii.gz' '"']; % output
cmd = [cmd ' -vector "' data_path filesep 'vect.nii.gz' '" -modulate none' ]; % output


cmd = [cmd ' "' data_path filesep 'DTIout.nii.gz' '"']; % output
end

[status, result] = msf_system(cmd);

%% Import reconstructed images metric
ADC = niftiread([data_path filesep 'adc.nii.gz']);
FA = niftiread([data_path filesep 'fa.nii.gz']);
VectorF = niftiread([data_path filesep 'vect.nii.gz']);
reD = VectorF(:,:,:,1);
greeN = VectorF(:,:,:,2);
bluE = VectorF(:,:,:,3);
DifT =0;

end

