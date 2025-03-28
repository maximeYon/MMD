function [status, result] =my_function_ANTs_Register(img_moving, img_out,Display)

%% Open NIFT dataset for MOCO
my_path_NIFTI = split(img_moving,filesep);
my_path_NIFTI = join(my_path_NIFTI(1:end-1),filesep);
my_path_NIFTI = my_path_NIFTI{1};

data_img_moving = squeeze(niftiread(img_moving));
infoNifti = niftiinfo(img_moving);

data_MOCO = zeros(size(data_img_moving));

%% Resave Nifti image for correction in loop
Nexp = size(data_img_moving,4);
moving_tmp_path = [my_path_NIFTI filesep 'moving_tmp.nii.gz'];
ref_tmp_path = [my_path_NIFTI filesep 'ref.nii.gz'];
out_tmp_path = [my_path_NIFTI filesep 'out_tmp.nii.gz'];

h = mdm_nii_h_empty;
h.pixdim(1:3) = infoNifti.PixelDimensions(1:3);
h.pixdim(4) = 1;
h.xyzt_units = 'SI';
h.dim_info = [1 2 3];
h.sform_code = 0; %

mdm_nii_write(data_img_moving(:,:,:,1), ref_tmp_path, h, 0);

wb = waitbar(0,'Motion correction in progress');
for ind = 1:Nexp
    mdm_nii_write(data_img_moving(:,:,:,ind), moving_tmp_path, h, 0);
 
    if ispc == 1 % Suppose a WSL installation
        cmd = 'bash -c -i "antsRegistration --verbose 1 --dimensionality 3';
        cmd = [cmd ' --output ['];
        cmd = [cmd '/mnt/' lower(moving_tmp_path(1)) strrep(moving_tmp_path(3:end),filesep,'/') ',']; % moving
        cmd = [cmd '/mnt/' lower(out_tmp_path(1)) strrep(out_tmp_path(3:end),filesep,'/') ']']; % out
        cmd = [cmd ' --interpolation Linear --winsorize-image-intensities [0.009,0.991]'];
        cmd = [cmd ' --initial-moving-transform ['];
        cmd = [cmd '/mnt/' lower(ref_tmp_path(1)) strrep(ref_tmp_path(3:end),filesep,'/') ',']; % ref is b0
        cmd = [cmd '/mnt/' lower(moving_tmp_path(1)) strrep(moving_tmp_path(3:end),filesep,'/') ',0]']; % moving
        cmd = [cmd ' --transform Rigid[0.08]'];
        cmd = [cmd ' --metric MI['];
        cmd = [cmd '/mnt/' lower(ref_tmp_path(1)) strrep(ref_tmp_path(3:end),filesep,'/') ',']; % ref is b0
        cmd = [cmd '/mnt/' lower(moving_tmp_path(1)) strrep(moving_tmp_path(3:end),filesep,'/') ',0.7,50,Regular,0.70]']; % moving
        cmd = [cmd ' --convergence [1000x500x5000x200,1e-6,30] --shrink-factors 8x4x2x1 --smoothing-sigmas 0.5x0x0x0vox --float 0 --collapse-output-transforms 1'];
    else
        %% to be programmed
    end
    
    if ispc == 1 % Suppose a WSL installation
        cmd = [cmd '"'];
        cmd = string(cmd);
    end
    
    [status, result] = msf_system(cmd);
    
    %     %% Alternative, copy to clipboard
%         reduced_cmd = cmd{1}(13:end-1);
%         clipboard('copy',reduced_cmd)
    
    
    data_MOCO(:,:,:,ind) = niftiread(out_tmp_path);
    waitbar((ind/Nexp),wb,'Motion correction in progress');
end
close(wb)

%% write corrected data;
h.sform_code = 0; %
h.pixdim(2:4) = infoNifti.PixelDimensions(1:3);
h.pixdim(1)=1;
mdm_nii_write(single(data_MOCO), img_out, h, 0);

%% Clean
delete(moving_tmp_path);
delete(ref_tmp_path);
delete(out_tmp_path);

if Display ==1
    %% Open Corrected and original nifti data
    data_noMOCO = niftiread(img_moving);
    % data_MOCO = niftiread(img_out);
    
    %% display Original and corrected images
    scliceN =round(size(data_noMOCO,3)/2);
    figure(4)
    for nrep = 1:min([size(data_noMOCO,4) 200])
        subplot(1,3,1)
        imagesc(squeeze(data_noMOCO(:,:,scliceN,nrep)))
        axis image
        title(['Original image n°' num2str(nrep)]);
        subplot(1,3,2)
        imagesc(squeeze(data_MOCO(:,:,scliceN,nrep)))
        axis image
        title(['Corrected image n°' num2str(nrep)]);
        subplot(1,3,3)
        imagesc(squeeze(data_noMOCO(:,:,scliceN,nrep))-squeeze(data_MOCO(:,:,scliceN,nrep)))
        axis image
        title('Difference');
        colormap gray
        drawnow
        pause(0.2)
    end
    
    sum_data_noMOCO = sum(data_noMOCO,4);
    sum_data_MOCO = sum(data_MOCO,4);
    figure(5)
    for nslice = 1:size(data_noMOCO,3)
        subplot(1,3,1)
        imagesc(squeeze(sum_data_noMOCO(:,:,nslice)))
        axis image
        title(['Original image, slice n°' num2str(nslice)]);
        subplot(1,3,2)
        imagesc(squeeze(sum_data_MOCO(:,:,nslice)))
        axis image
        title(['Corrected image, slice n°' num2str(nslice)]);
        subplot(1,3,3)
        imagesc(squeeze(sum_data_noMOCO(:,:,nslice))-squeeze(sum_data_MOCO(:,:,nslice)))
        axis image
        title('Difference');
        colormap gray
        drawnow
        pause(2)
    end
    
end
