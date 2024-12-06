%Run bootstrap analysis
clear all

wd = pwd;

%method = 'dti_euler';
method = 'dtd';
%method = 'dtr1d';

%Define paths to bootstrap folders
% datasets_paths{1} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd08/mdd';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd07/mdd';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd06/mdd';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd05/mdd';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd04/mdd';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd03/mdd';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd02/mdd';
% % datasets_paths{1} = '/Users/daniel/Dropbox/NMRdata/Caeyenberghs/FWF_DTC01_Caeyenberghs/mdd';
% % datasets_paths{1} = '/Users/daniel/Dropbox/NMRdata/United/data4dtdpaper';
% % datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191022/mdd';
% % datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191023/mdd';
% % datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191023_pm/mdd';
% % datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191029pm_v1/mdd';
% % datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191029night_v1/mdd';
% % datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191030am_p1/mdd';
% % datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd';
% % datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191109_p04/mdd';
% 
% Ndata = numel(datasets_paths);
% in_paths = cell(0,0);
% bs_paths = cell(0,0);
% out_paths = cell(0,0);
% for ndata = 1:Ndata
%     datasets_path = datasets_paths{ndata};
%     expnams = mdm_bruker_dir2expnams(datasets_path);
%     Nexp = numel(expnams);
%     for nexp = 1:Nexp
%         expnam = expnams{nexp};
%         in_paths{1+numel(in_paths)} = fullfile(datasets_path,expnam,'nii_xps');
%         bs_paths{1+numel(bs_paths)} = fullfile(datasets_path,expnam,method,'bootstraps');
%         out_paths{1+numel(out_paths)} = fullfile(datasets_path,expnam,method,'maps');
%     end
% end
in_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/nii_xps'};
bs_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/dtd/bootstraps'};
out_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/dtd/maps'};

% clear in_paths bs_paths
% method = 'dtr2d';
% in_paths{1} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_15_20190725_161749/mdd/nii_xps_merge_21to30';
% in_paths{1+numel(in_paths)} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_17_20190726_115920/mdd/nii_xps_merge_21to30';
% in_paths{1+numel(in_paths)} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_14_20190725_122035/mdd/nii_xps_merge';
% Ndata = numel(in_paths);
% for ndata = 1:Ndata
%     bs_paths{ndata} = fullfile(in_paths{ndata},method,'bootstraps');
% end


nii_name = 'data';


Ndata = numel(in_paths);
for ndata = 1:Ndata
    in_path = in_paths{ndata};
    bs_path = bs_paths{ndata};
    out_path = out_paths{ndata};
    msf_mkdir(out_path);

    % Connect to data
    clear s
    s.nii_fn = fullfile(in_path, [nii_name '.nii.gz']);
    s.mask_fn = fullfile(in_path, [nii_name '_mask.nii.gz']);
    s.xps = mdm_xps_load(fullfile(in_path, [nii_name '_xps.mat']));

    [I_in,h]   = mdm_nii_read(s.nii_fn);
    M = mdm_mask_load(s);
    xps = s.xps;

    I_in = double(I_in).*double(M);
    I_in = permute(I_in,[2 1 3 4]);
    I_in = flip(I_in,1);
    sz_tile = size(I_in);
    pixaspect = h.pixdim(3)/h.pixdim(2);
    imaspect = sz_tile(2)/sz_tile(1);
    
    nk = round(sz_tile(3)/2);

    Imax = max(I_in(:));
    im2d_in = reshape(I_in(:,:,nk,:),[sz_tile(1) sz_tile(2) 1 sz_tile(4)]);
    
    bsno = msf_getdirno(bs_path);        
    bs_dps = cell(numel(bsno),1);
    
    im2d_dif_array = zeros(sz_tile(1),sz_tile(2),1,sz_tile(4),numel(bsno));
    
    parfor nbs = 1:numel(bsno)
        o     = fullfile(bs_path,num2str(nbs));
        mfs_fn   = fullfile(o, 'mfs.mat');
        if exist(mfs_fn,'file')==2

            mfs = mdm_mfs_load(mfs_fn);

            m = mfs.m;        
            sz = size(m);

            I_fit = zeros(sz(1),sz(2),sz(3),xps.n);
            %for nk = 1:sz(3)
                for nj = 1:sz(2)
                    for ni = 1:sz(1)
                        if M(ni,nj,nk)
                            %s_fit = dtr2d_1d_fit2data(squeeze(m(ni,nj,nk,:))', xps);
                            s_fit = dtd_1d_fit2data(squeeze(m(ni,nj,nk,:))', xps);
                            I_fit(ni,nj,nk,:) = s_fit;
                            %figure(1), clf, plot(1:xps.n,squeeze(I_fit(ni,nj,nk,:)),'.'), pause(.1)
                        end
                    end
                end
            %end

            I_fit = double(I_fit).*double(M);
            I_fit = permute(I_fit,[2 1 3 4]);
            I_fit = flip(I_fit,1);

            I_dif = I_in - I_fit;

            im2d_dif = reshape(I_dif(:,:,nk,:),[sz_tile(1) sz_tile(2) 1 sz_tile(4)]);

            im2d_dif_array(:,:,:,:,nbs) = im2d_dif;
        end
    end
    
    im2d_dif = sum(im2d_dif_array,5)/numel(bsno);
    
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

    print(fullfile(out_path,'qualitycontrol'),'-loose','-dpdf')

end

