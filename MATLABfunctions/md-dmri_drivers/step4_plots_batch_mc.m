%Plot parameters maps and global stats
clear all

wd = pwd;

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

% in_paths = {'/Users/daniel/Dropbox/NMRdata/md-dmri-data/topgaard_nmrbiomed2019/SmallSpheres_BigSpheres_OrderedSticks/pdata_mddmri/indata'};
% bs_paths = {'/Users/daniel/Dropbox/NMRdata/md-dmri-data/topgaard_nmrbiomed2019/SmallSpheres_BigSpheres_OrderedSticks/pdata_mddmri/fitdata/dtd/bootstrap'};
% out_paths = {'/Users/daniel/Dropbox/NMRdata/md-dmri-data/topgaard_nmrbiomed2019/SmallSpheres_BigSpheres_OrderedSticks/pdata_mddmri/fitdata/maps'};
% in_paths = {'/Users/daniel/Dropbox/NMRdata/md-dmri-data/topgaard_nmrbiomed2019/SmallSpheres_BigSpheres_RandomSticks/pdata_mddmri/indata'};
% bs_paths = {'/Users/daniel/Dropbox/NMRdata/md-dmri-data/topgaard_nmrbiomed2019/SmallSpheres_BigSpheres_RandomSticks/pdata_mddmri/fitdata/dtd/bootstrap'};
% out_paths = {'/Users/daniel/Dropbox/NMRdata/md-dmri-data/topgaard_nmrbiomed2019/SmallSpheres_BigSpheres_RandomSticks/pdata_mddmri/fitdata/maps'};
% in_paths = {'/Users/daniel/Dropbox/NMRdata/Cardiff/HongPhantom/mdd/MDD_372/nii_xps'};
% bs_paths = {'/Users/daniel/Dropbox/NMRdata/Cardiff/HongPhantom/mdd/MDD_372/pdata/dtd/bootstrap'};
% out_paths = {'/Users/daniel/Dropbox/NMRdata/Cardiff/HongPhantom/mdd/MDD_372/pdata/dtd/maps'};

% in_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/MR750_20190521/RWI_volunteer/mdd/MDD_intermediate/pdata/dtd/bootstrap'};
% bs_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/MR750_20190521/RWI_volunteer/mdd/MDD_intermediate/pdata/dtd/bootstrap'};
% out_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/MR750_20190521/RWI_volunteer/mdd/MDD_intermediate/pdata/pdf_maps'};

% in_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/MR750_20190521/RWI_volunteer/mdd/MDD_intermediate/pdata/dtd/bootstrap'};
% bs_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/MR750_20190521/RWI_volunteer/mdd/MDD_intermediate/pdata/dtd/bootstrap'};
% out_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/MR750_20190521/RWI_volunteer/mdd/MDD_intermediate/pdata/pdf_maps'};
in_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/nii_xps'};
bs_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/dtd/bootstraps'};
out_paths = {'/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/mc/dtd/maps'};


% 
% clear in_paths bs_paths
% method = 'dtr2d';
% in_paths{1} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_15_20190725_161749/mdd/nii_xps_merge_21to30';
% in_paths{1+numel(in_paths)} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_17_20190726_115920/mdd/nii_xps_merge_21to30';
% in_paths{1+numel(in_paths)} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_14_20190725_122035/mdd/nii_xps_merge';
% Ndata = numel(in_paths);
% for ndata = 1:Ndata
%     bs_paths{ndata} = fullfile(in_paths{ndata},method,'bootstraps');
%     out_paths{ndata} = fullfile(bs_paths{ndata},'maps');
% end
% 

% clear in_paths bs_paths
% method = 'dtr2d';
% bs_paths = {'/Users/daniel/Dropbox/NMRdata/Cardiff/Prisma/June2018/bootstraps_mc_eddyTopup'};
% out_paths = {'/Users/daniel/Dropbox/NMRdata/Cardiff/Prisma/June2018/pdfmaps_mc_eddyTopup'};


%Define bins
disomin = [0 0 2.5]*1e-9; disomax = [2.5 2.5 5]*1e-9;
dratiomin = [1 1 1]*eps; dratiomax = [1 1 1]/eps;
sddeltamin = [.25 0 0]; sddeltamax = [1 .25 1];
% disomin = [0 0 0]*1e-9; disomax = [5 5 5]*1e-9;
% dratiomin = [4 eps 1/4]; dratiomax = [1/eps 1/4 4];
% sddeltamin = [0 0 0]; sddeltamax = [1 1 1];
% disomin = [0 0 .5]*1e-9; disomax = [.5 .5 1.5]*1e-9;
% dratiomin = [1 1 1]*eps; dratiomax = [1 1 1]/eps;
% sddeltamin = [.25 0 0]; sddeltamax = [1 .25 1];
% disomin = [0 0 1 2.5]*1e-9; disomax = [2.5 1 2.5 5]*1e-9;
% dratiomin = [1 1 1 1]*eps; dratiomax = [1 1 1 1]/eps;
% sddeltamin = [.25 0 0 0]; sddeltamax = [1 .25 .25 1];
%r2min = [1 1 1]; r2max = 80*[1 1 1];
r2min = .5*[1 1 1]; r2max = 30*[1 1 1];
r1min = .1*[1 1 1]; r1max = 2*[1 1 1];

%Define color limits for parameter maps
clim.s0 = .5*[0 1]; % Multiplied with s0 below
clim.mdiso = 3.5e-9*[0 1]; %clim.mdiso = 1e-9*[0 1];
clim.msddelta = 1*[0 1];
clim.mr2 = max(r2max)*[0 1];
clim.mr1 = max(r1max)*[0 1]; clim.mr1 = .8*[0 1];
clim.vdiso = .3*3e-9^2*[0 1]; %clim.vdiso = .3*1e-9^2*[0 1];
clim.vsddelta = .2*[0 1];
clim.vr2 = .2*max(r2max)^2*[0 1];
clim.vr1 = .2*max(r1max)^2*[0 1]; clim.vr1 = .2*.8^2*[0 1];
clim.cvdisosddelta = .1*3e-9*1*[-1 1]; %clim.cvdisosddelta = .2*1e-9*1*[-1 1];
clim.cvdisor2 = .1*max(r2max)*3e-9*[-1 1];
clim.cvsddeltar2 = .1*max(r2max)*1*[-1 1];
clim.cvdisor1 = .1*max(r1max)*3e-9*[-1 1]; clim.cvdisor1 = .1*.8*3e-9*[-1 1];
clim.cvsddeltar1 = .1*max(r1max)*1*[-1 1]; clim.cvsddeltar1 = .1*.8*1*[-1 1];
clim.mask_threshold = .01;

%------------------------------

opt.(method).bin_disomin = disomin; opt.(method).bin_disomax = disomax;
opt.(method).bin_dratiomin = dratiomin; opt.(method).bin_dratiomax = dratiomax;
opt.(method).bin_sddeltamin = sddeltamin; opt.(method).bin_sddeltamax = sddeltamax;
if strcmp(method,'dtr2d')
    opt.(method).bin_r2min = r2min; opt.(method).bin_r2max = r2max;
elseif strcmp(method,'dtr1d')
    opt.(method).bin_r1min = r1min; opt.(method).bin_r1max = r1max;
end

Ndata = numel(in_paths);
tic
for ndata = 1:Ndata

    bs_path = bs_paths{ndata};
    out_path = out_paths{ndata};

    bs_dps = mdm_dps_collectbs(method, bs_path, opt);

        
    if ~all(cellfun('isempty',bs_dps))
        median_dps = mdm_dps_median(bs_dps);
        %%
        clim.mask_threshold = .03;
        mplot_technicolor(method, median_dps, fullfile(out_path,'technicolor'), clim)
        mplot_technicolor_slices(method, median_dps, fullfile(out_path,'slices'), clim)
        mplot_technicolor_allmapsperslice(method, median_dps, fullfile(out_path,'slices'), clim)
        mplot_globalstats(method, median_dps, fullfile(out_path,'globalstats'), clim)
%%
    end
   
end
toc


