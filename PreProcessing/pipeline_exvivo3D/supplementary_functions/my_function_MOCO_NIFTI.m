function my_function_MOCO_NIFTI(pathData,mode)
%% Open topUp dataset for MOCO
img = squeeze(niftiread(pathData));

%% parameter: slice reference number
Nslice = round(size(img,3)/2); % slice reference number

%% Get reference image by averaging the first 5 images
ref_img_up = sum(img(:,:,:,1:5),4);

%% display ref images
% figure(1)
% imagesc(ref_img(:,:,Nslice)')
% colormap gray

%% Calculate translation of the image serie in a 2D slice
[optimizer,metric] = imregconfig("multimodal");
displacementUp = zeros(2,size(img,4));

for ind_img = 1:size(img,4)
    % [moving_reg,R_reg] = imregister(squeeze(img_up(:,:,Nslice,ind_img)),ref_img_up(:,:,Nslice),"translation",optimizer,metric); % "rigid" to add rotation
    tform = imregtform(squeeze(img(:,:,Nslice,ind_img)),ref_img_up(:,:,Nslice),"translation",optimizer,metric);
    displacement(:,ind_img) = tform.T(3,1:2);
end
% figure(2)
% subplot(1,2,1)
% imshowpair(ref_img_up(:,:,Nslice)',squeeze(img_up(:,:,Nslice,ind_img))',"Scaling","joint")
% title('Without registration')
% subplot(1,2,2)
% imshowpair(ref_img_up(:,:,Nslice)',moving_reg',"Scaling","joint")
% title('With registration')

%% remove outsider from displacement curves
indOut = zeros(2,size(img,4));
for ind_disp = 1:2
    [~,indOut(ind_disp,:)] = rmoutliers(displacement(ind_disp,:),"gesd");
end
indOut = logical(sum(indOut,1));


%% filter displacement curves
x = 1:size(img,4);
w = abs(indOut-1);

displacement_filtered = csaps(x,displacement,0.000001,x,w);
Mean_displacement = displacement_filtered';

%% display
figure(3)
hold on
plot(displacement','r')
plot(displacement_filtered','k')
plot(indOut'*std(sum(displacement_filtered,1))*30-std(sum(displacement_filtered,1))*10,'b')
% ylim([min(indOutUp'*std(sum(displacementUp_filtered,1))*30-std(sum(displacementUp_filtered,1))*12) max(indOutUp'*std(sum(displacementUp_filtered,1))*30-std(sum(displacementUp_filtered,1))*8)])
legend('raw', '','filtered','','outliers')
title('motion curves')
ylabel('motion in pixels');
drawnow

data_path = split(pathData,filesep);
data_path = join(data_path(1:end-1,1),filesep,1); data_path = data_path{1};
save([data_path filesep 'MotionCurves.mat'],'displacement','displacement_filtered','indOut');

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
   if max(max(abs(Mean_displacement))) > 0.7
       answer = 'Yes';
       saveas(gcf,[data_path filesep 'Motion_curves_' answer 'MOCO.jpeg'])
   else
       answer = 'No';
       saveas(gcf,[data_path filesep 'Motion_curves_' answer 'MOCO.jpeg'])
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
        img_size =size(squeeze(img(:,:,1,1)));
        img_corr = img;
        for ind_img = 1:size(img,4)
            my_T = [1 0 0; 0 1 0;0 0 1];
            my_T(3,1:2) = Mean_displacement(ind_img,:);
            tform = affine2d(my_T);
            for ind_slice = 1:size(img,3)
%                 img_corr(:,:,ind_slice,ind_img) = imwarp(img(:,:,ind_slice,ind_img),tform,"OutputView",imref2d(img_size),"interp","linear");
                img_corr(:,:,ind_slice,ind_img) = abs(fineshift(img(:,:,ind_slice,ind_img),-tform.T(3,[2 1])));
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
        my_save_NIFTI(img_corr,data_path_pv,[data_path filesep 'dataMOCO.nii.gz']);
        save([data_path filesep 'MOCO_UpDown_displacement.mat'],'Mean_displacement');
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




