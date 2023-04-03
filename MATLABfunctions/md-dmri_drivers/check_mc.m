%Run bootstrap analysis
clear all

wd = pwd;

% in_fn = '/Users/daniel/Dropbox/NMRdata/Caeyenberghs/FWF_DTC01_Caeyenberghs/mdd/MDD_82_nob0/nii_xps/data.nii.gz';
% mc_fn = '/Users/daniel/Dropbox/NMRdata/Caeyenberghs/FWF_DTC01_Caeyenberghs/mdd/MDD_82_fsl_nob0/nii_xps/data.nii.gz';
% out_fn = '/Users/daniel/Dropbox/NMRdata/Caeyenberghs/FWF_DTC01_Caeyenberghs/mdd/MDD_82_fsl_nob0/mc_check.pdf';

in_fn = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/mc_check/MDD_intermediate_FS.nii.gz';
mc_fn = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/mc_check/MDD_intermediate_FS_mc.nii.gz';
out_fn = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/mc_check//mc_check.pdf';

[I_in,h]   = mdm_nii_read(in_fn);
[I_mc,~]   = mdm_nii_read(mc_fn);

I_in = double(I_in);
I_in = permute(I_in,[2 1 3 4]);
I_in = flip(I_in,1);
sz_tile = size(I_in);
pixaspect = h.pixdim(3)/h.pixdim(2);
imaspect = sz_tile(2)/sz_tile(1);
    
nk = round(sz_tile(3)/2);

Imax = max(I_in(:));
im2d_in = reshape(I_in(:,:,nk,:),[sz_tile(1) sz_tile(2) 1 sz_tile(4)]);
    
I_mc = double(I_mc);
I_mc = permute(I_mc,[2 1 3 4]);
I_mc = flip(I_mc,1);

im2d_mc = reshape(I_mc(:,:,nk,:),[sz_tile(1) sz_tile(2) 1 sz_tile(4)]);

im2d_dif = im2d_in - im2d_mc;

map_in = colormap(gray(64));
Ncindex_in = size(map_in,1);
map_dif = mplot_cmaphotcold(64);
Ncindex_dif = size(map_dif,1);

im2d_in_indexed = uint8(Ncindex_in*(10*im2d_in/Imax));
im2d_dif_indexed = uint8(Ncindex_dif*(10*im2d_dif/2/Imax+.5));
out_in = imtile(im2d_in_indexed,map_in);
out_dif = imtile(im2d_dif_indexed,map_dif);
figure(1), clf
axes('position',[0 0 .5 1])
imshow(out_in);
axes('position',[.5 0 .5 1])
imshow(out_dif);

papersize = 15*[2 1];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

print(out_fn,'-loose','-dpdf')



