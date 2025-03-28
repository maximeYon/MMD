function my_function_MOCO_NIFTI_rigid(pathData,mode)
%% Open topUp dataset for MOCO
img = squeeze(niftiread(pathData));

%% parameter: slice reference number
Nslice = round(size(img,3)/2); % slice reference number

%% Get reference image by averaging the first 5 images
% ref_img_up = sum(img(:,:,:,1:5),4);
ref_img_up = sum(img./repmat(reshape(mean(reshape(img,size(img,1)*size(img,2)*size(img,3),size(img,4)),1),1,1,1,size(img,4)),size(img,1),size(img,2),size(img,3),1),4);
mask = zeros(size(ref_img_up));
mask(ref_img_up>mean(ref_img_up(:))./0.75)=1;
%% display ref images
% % figure(1)
% subplot(1,2,1)
% imagesc(ref_img_up(:,:,Nslice)')
% colormap gray
% subplot(1,2,2)
% imagesc(mask(:,:,Nslice)')
% colormap gray

%% Calculate translation of the image serie in a 2D slice
[optimizer,metric] = imregconfig("multimodal");
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;
translation = zeros(1,2,size(img,4));
rotation = zeros(1,size(img,4));

wb = waitbar(0,'Motion correction in progress');
for ind_img = 1:size(img,4)
    % [moving_reg,R_reg] = imregister(squeeze(img_up(:,:,Nslice,ind_img)),ref_img_up(:,:,Nslice),"translation",optimizer,metric); % "rigid" to add rotation
    tform = imregtform(squeeze(img(:,:,Nslice,ind_img)).*squeeze(mask(:,:,Nslice)),ref_img_up(:,:,Nslice).*squeeze(mask(:,:,Nslice)),"rigid",optimizer,metric);
    structT{ind_img} = tform;
    translation(:,:,ind_img) = tform.Translation;
    rotation(:,ind_img) = tform.RotationAngle;
    waitbar((ind_img/size(img,4)),wb,'Motion correction in progress');
end
close(wb)
translation = squeeze(translation);

%% remove outsider from displacement curves
indOut = zeros(3,size(img,4));
for ind_disp = 1:2
    [~,indOut(ind_disp,:)] = rmoutliers(translation(ind_disp,:),"gesd");
    [~,indOut(3,:)] = rmoutliers(rotation,"gesd");
end
indOut = logical(sum(indOut,1));

figure(1)
subplot(3,1,1)
plot(rotation)
title('Rotation')
subplot(3,1,2)
plot(translation')
title('Translation')
subplot(3,1,3)
plot(indOut)
title('Outliers')

%% Ask user for correction or not
if strcmp(mode,'auto')
    if sum(indOut)<1
        answer = 'Yes';
        saveas(gcf,[data_path filesep 'Motion_curves_' answer 'MOCO.jpeg'])
    else
        answer = 'No';
        saveas(gcf,[data_path filesep 'Motion_curves_' answer 'MOCO.jpeg'])
    end
else
    answer = questdlg('Do you want to apply Motion correction?', ...
        'Apply Motion Correction', ...
        'Yes','No','Yes');
    % Handle response
end
switch answer
    case 'Yes'
        disp([answer 'Applying MOCO'])
        %% Apply displacement to the image serie
        img_size =size(squeeze(img(:,:,1,1)));
        img_corr = img;
        for ind_img = 1:size(img,4)
            tform = structT{ind_img};
            for ind_slice = 1:size(img,3)
                img_corr(:,:,ind_slice,ind_img) = imwarp(img(:,:,ind_slice,ind_img),tform,"OutputView",imref2d(img_size),"interp","linear");
                % img_corr(:,:,ind_slice,ind_img) = abs(fineshift(img(:,:,ind_slice,ind_img),-tform.T(3,[2 1])));
            end
        end

        aa =1;
        figure(2)
        for ind_displ = 1:size(img,4)
            subplot(1,2,1)
            % imshowpair(ref_img_up(:,:,Nslice)',squeeze(img(:,:,Nslice,ind_displ))'.*10,"Scaling","joint")
            imagesc(squeeze(img(:,:,Nslice,ind_displ))')
            colormap gray
            title(['Without registration n° ' num2str(ind_displ)])
            subplot(1,2,2)
            % imshowpair(ref_img_up(:,:,Nslice)',squeeze(img_corr(:,:,Nslice,ind_displ))'.*10,"Scaling","joint")
            imagesc(squeeze(img_corr(:,:,Nslice,ind_displ))')
            colormap gray
            title(['With registration n° ' num2str(ind_displ)])
            drawnow
            pause(0.2)
        end


        img_sum = img./repmat(reshape(mean(reshape(img,size(img,1)*size(img,2)*size(img,3),size(img,4)),1),1,1,1,size(img,4)),size(img,1),size(img,2),size(img,3),1);
        img_corr_sum = img_corr./repmat(reshape(mean(reshape(img_corr,size(img,1)*size(img,2)*size(img,3),size(img,4)),1),1,1,1,size(img,4)),size(img,1),size(img,2),size(img,3),1);

         figure(3)
          subplot(1,2,1)
           imagesc(sum(img_sum(:,:,Nslice,:),4)')
            colormap gray
            title('Without registration')
             subplot(1,2,2)
           imagesc(sum(img_corr_sum(:,:,Nslice,:),4)')
            colormap gray
            title('Without registration')

        %% Save corrected NIFTI files Up and Down
        data_path_pv = split(pathData,filesep);
        data_path_pv = join(data_path_pv(1:end-2,1),filesep,1); data_path_pv = data_path_pv{1};
        my_save_NIFTI(img_corr,data_path_pv,[data_path filesep 'dataMOCO.nii.gz']);
        save([data_path filesep 'MOCO_UpDown_displacement.mat'],'rotation','translation','indOut');
        disp('data MOCO saved');

    case 'No'
        disp([answer 'Bye'])
end

end

%% function save NIFTI
function my_save_NIFTI(data,data_path_pv,nii_fn)
PVM_SpatResol=ReadPV360Param([data_path_pv filesep], 'PVM_SpatResol') ;
PVM_SliceThick=ReadPV360Param([data_path_pv filesep], 'PVM_SliceThick') ;
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
h.dim_info = [1 2 3];
h.sform_code = 0;
mdm_nii_write(data, nii_fn, h, 0);
end




