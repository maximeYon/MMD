function my_function_MOCO_NIFTI_topup(pathData,mode)
%% Open topUp dataset for MOCO

%% Denoised NIFTI: Import reconstructed imagesmode
if isfile([pathData filesep 'dataUp.nii.gz'])
    pathDataUp = [pathData filesep 'dataUp.nii.gz'];
    pathDataDown = [pathData filesep 'dataDown.nii.gz'];
    in_denoise = 0;
else
    pathDataUp = [pathData filesep 'denoised' filesep 'dataUp.nii.gz'];
    pathDataDown = [pathData filesep 'denoised' filesep 'dataDown.nii.gz'];
    in_denoise = 1;
end
img_up = squeeze(niftiread(pathDataUp));
img_down = squeeze(niftiread(pathDataDown));

%% parameter: slice reference number
Nslice = round(size(img_up,3)/2); % slice reference number

%% Get reference image by averaging the first 5 images
ref_img_up = sum(img_up(:,:,:,1:5),4);
ref_img_down = sum(img_down(:,:,:,1:5),4);

%% display ref images
% figure(1)
% subplot(1,2,1)
% imagesc(ref_img_up(:,:,Nslice)')
% subplot(1,2,2)
% imagesc(ref_img_down(:,:,Nslice)')
% colormap gray

%% Calculate translation of the image serie in a 2D slice
[optimizer,metric] = imregconfig("multimodal");
displacementUp = zeros(2,size(img_up,4));
displacementDown = zeros(2,size(img_up,4));
for ind_img = 1:size(img_up,4)
    % [moving_reg,R_reg] = imregister(squeeze(img_up(:,:,Nslice,ind_img)),ref_img_up(:,:,Nslice),"translation",optimizer,metric); % "rigid" to add rotation
    tformUp = imregtform(squeeze(img_up(:,:,Nslice,ind_img)),ref_img_up(:,:,Nslice),"translation",optimizer,metric);
    displacementUp(:,ind_img) = tformUp.T(3,1:2);
    tformDown = imregtform(squeeze(img_down(:,:,Nslice,ind_img)),ref_img_down(:,:,Nslice),"translation",optimizer,metric);
    displacementDown(:,ind_img) = tformDown.T(3,1:2);
end
% figure(2)
% subplot(1,2,1)
% imshowpair(ref_img_up(:,:,Nslice)',squeeze(img_up(:,:,Nslice,ind_img))',"Scaling","joint")
% title('Without registration')
% subplot(1,2,2)
% imshowpair(ref_img_up(:,:,Nslice)',moving_reg',"Scaling","joint")
% title('With registration')

%% remove outsider from displacement curves
indOutUp = zeros(2,size(img_up,4));
indOutDown = zeros(2,size(img_up,4));
for ind_disp = 1:2
    [~,indOutUp(ind_disp,:)] = rmoutliers(displacementUp(ind_disp,:),"gesd");
    [~,indOutDown(ind_disp,:)] = rmoutliers(displacementDown(ind_disp,:),"gesd");
end
indOutUp = logical(sum(indOutUp,1));
indOutDown = logical(sum(indOutDown,1));

%% filter displacement curves
x = 1:size(img_up,4);
wUp = abs(indOutUp-1);
wDown = abs(indOutUp-1);
displacementUp_filtered = csaps(x,displacementUp,0.000001,x,wUp);
displacementDown_filtered = csaps(x,displacementDown,0.000001,x,wDown);
Mean_displacement = (displacementUp_filtered' + displacementDown_filtered')/2;

%% display
figure(3)
subplot(1,3,1)
hold on
plot(displacementUp','r')
plot(displacementUp_filtered','k')
plot(indOutUp'*std(sum(displacementUp_filtered,1))*30-std(sum(displacementUp_filtered,1))*10,'b')
% ylim([min(indOutUp'*std(sum(displacementUp_filtered,1))*30-std(sum(displacementUp_filtered,1))*12) max(indOutUp'*std(sum(displacementUp_filtered,1))*30-std(sum(displacementUp_filtered,1))*8)])
legend('raw', '','filtered','','outliers')
title('motion blip up')
ylabel('motion in pixels');
subplot(1,3,2)
hold on
plot(displacementDown','r')
plot(displacementDown_filtered','k')
plot(indOutDown'*std(sum(displacementDown_filtered,1))*30-std(sum(displacementDown_filtered,1))*10,'b')
% ylim([min(indOutUp'*std(sum(displacementDown_filtered,1))*30-std(sum(displacementDown_filtered,1))*12) max(indOutUp'*std(sum(displacementDown_filtered,1))*30-std(sum(displacementDown_filtered,1))*8)])
legend('raw', '','filtered','','outliers')
title('motion blip Down')
ylabel('motion in pixels');
subplot(1,3,3)
hold on
plot(displacementUp_filtered','r')
plot(displacementDown_filtered','b')
plot(Mean_displacement,'k')
legend('diplacement Up filtered', '','diplacement Down filtered','','Mean displacement filtered','Location','best')
drawnow

diff_disp = abs(displacementUp_filtered-displacementDown_filtered);

save([pathData filesep 'MotionCurves.mat'],'displacementUp','displacementUp_filtered','displacementDown','displacementDown_filtered','indOutUp','indOutDown');

% EpiEffBandwidth=ReadPV360Param([pathData filesep], 'PVM_EpiEffBandwidth') ;
% EffSWh=ReadPV360Param([pathData filesep], 'PVM_EffSWh') ;
% Matrix=ReadPV360Param([pathData filesep], 'PVM_Matrix') ;
% figure(4)
% hold on
% plot((Mean_displacement(:,1)/Matrix(1,1))*EffSWh,'b')
% plot((Mean_displacement(:,2)/Matrix(1,2))*EpiEffBandwidth,'r')
% legend('drift in read', 'drift in phase');
% ylabel('drift in Hz');

%% Ask user for correction or not
if strcmp(mode,'auto')
   if max(max(abs(Mean_displacement))) > 0.7 && max(max(abs(diff_disp)))< max(max(abs(Mean_displacement)))/4
       answer = 'Yes';
       saveas(gcf,[pathData filesep 'Motion_curves_' answer 'MOCO.jpeg'])
   else
       answer = 'No';
       saveas(gcf,[pathData filesep 'Motion_curves_' answer 'MOCO.jpeg'])
   end
else
answer = questdlg('Do you want to apply Motion correction?', ...
    'Aply Motion Correction', ...
    'Yes','No','Yes');
% Handle response
end
switch answer
    case 'Yes'
        disp([answer 'Applying MOCO'])
        %% Apply displacement to the image serie
        img_size =size(squeeze(img_up(:,:,1,1)));
        img_up_corr = img_up;
        img_down_corr = img_down;
        for ind_img = 1:size(img_up,4)
            my_T = [1 0 0; 0 1 0;0 0 1];
            my_T(3,1:2) = Mean_displacement(ind_img,:);
            tform = affine2d(my_T);
            for ind_slice = 1:size(img_up,3)
%                 img_up_corr(:,:,ind_slice,ind_img) = imwarp(img_up(:,:,ind_slice,ind_img),tform,"OutputView",imref2d(img_size),"interp","linear");
%                 img_down_corr(:,:,ind_slice,ind_img) = imwarp(img_down(:,:,ind_slice,ind_img),tform,"OutputView",imref2d(img_size),"interp","linear");                
                img_up_corr(:,:,ind_slice,ind_img) = abs(fineshift(img_up(:,:,ind_slice,ind_img),-tform.T(3,[2 1])));
                img_down_corr(:,:,ind_slice,ind_img) = abs(fineshift(img_down(:,:,ind_slice,ind_img),-tform.T(3,[2 1])));
            end
        end
%         
%         %% Display
%         figure(11)
%         clims = [0 max(max((sum(img_up(:,:,3,:),4))))];
%         subplot(3,3,1)
% %         imagesc(flip(img_up(:,:,3,129)',1))
%         imagesc(flip(sum(img_up(:,:,3,:),4)',1),clims)
%         axis image
%         
%         subplot(3,3,4)
% %         imagesc(flip(img_up_corr(:,:,3,129)',1))
%         imagesc(flip(sum(img_up_corr(:,:,3,:),4)',1),clims)
%         axis image
%         subplot(3,3,5)
% %         imagesc(flip(img_up_corr2(:,:,3,129)',1))
%         imagesc(flip(sum(img_up_corr2(:,:,3,:),4)',1),clims)
%         axis image
%         subplot(3,3,6)
% %         imagesc(flip(img_up_corr(:,:,3,129)',1)-flip(img_up_corr2(:,:,3,129)',1))
%         imagesc(flip(sum(img_up_corr(:,:,3,:),4)',1)-flip(sum(img_up_corr2(:,:,3,:),4)',1))
%         colorbar
%         colormap gray
%         axis image
%         
%         subplot(3,3,7)
% %         imagesc(flip(img_up(:,:,3,129)',1)-flip(img_up_corr(:,:,3,129)',1))
%         imagesc(flip(sum(img_up(:,:,3,:),4)',1)-flip(sum(img_up_corr(:,:,3,:),4)',1))
%         colorbar
%         colormap gray
%         axis image
%         subplot(3,3,8)
% %         imagesc(flip(img_up(:,:,3,129)',1)-flip(img_up_corr2(:,:,3,129)',1))
%         imagesc(flip(sum(img_up(:,:,3,:),4)',1)-flip(sum(img_up_corr2(:,:,3,:),4)',1))
%         colorbar
%         colormap gray
%         axis image
%         linkaxes
        
        
        %         figure(5)
        %         ind_displ =390;
        %         subplot(1,2,1)
        %         imshowpair(ref_img_up(:,:,Nslice)',squeeze(img_up(:,:,Nslice,ind_displ))',"Scaling","joint")
        %         title('Without registration')
        %         subplot(1,2,2)
        %         imshowpair(ref_img_up(:,:,Nslice)',squeeze(img_up_corr(:,:,Nslice,ind_displ))',"Scaling","joint")
        %         title('With registration')
        
        %% Save corrected NIFTI files Up and Down
        data_path_pv = split(pathData,filesep);
        data_path_pv = join(data_path_pv(1:end-2,1),filesep,1); data_path_pv = data_path_pv{1};
        my_save_NIFTI(img_up_corr,data_path_pv,[pathData filesep 'dataUp.nii.gz']);
        my_save_NIFTI(img_down_corr,data_path_pv,[pathData filesep 'dataDown.nii.gz']);
        save([pathData filesep 'MOCO_UpDown_displacement.mat'],'Mean_displacement');
        disp('dataUp and dataDown saved');
        
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




