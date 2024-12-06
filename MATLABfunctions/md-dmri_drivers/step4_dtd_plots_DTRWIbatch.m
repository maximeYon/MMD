%Plot parameters maps and global stats
clear all

% Define full paths to the nifti files to be analyzed 
% data_path =  '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419';
data_path =  '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
% data_dir = dir(fullfile(data_path,'*'));
% data_dir = dir(fullfile(data_path,'P025237*')); opt.k_range = 12;
data_dir = dir(fullfile(data_path,'P0240306*')); opt.k_range = 15;
% data_dir = dir(fullfile(data_path,'P0275752*')); opt.k_range = 13;
% data_dir = dir(fullfile(data_path,'P0333068*')); opt.k_range = 16;
% data_dir = dir(fullfile(data_path,'P0297359*')); opt.k_range = 16;
% data_dir = dir(fullfile(data_path,'P0275680*')); opt.k_range = 9;

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
% % data_path = 'C:/Users/topga/Dropbox/NMRdata/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'2020*'));

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_Mahajan_1/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'20*'));
% 
% data_path = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'20*'));
% 
% data_path = '/Users/daniel/Dropbox/NMRdata/Philips_QC/acquisition_date_time';
% data_path = '/Users/daniel/Dropbox/NMRdata/GE_QC/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'*'));

% data_path = '/Users/daniel/Dropbox/NMRdata/Philips_Sahlgrenska/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'*'));

nii_fns = cell(size(data_dir,1),1);
for ndata = 1:size(data_dir,1)
    nii_fns{ndata} = fullfile(data_path,data_dir(ndata).name,'rwi','nii_xps','data_mc.nii.gz');
end
% nii_fns = nii_fns(end:-1:1);
% nii_fns = nii_fns(71);
% nii_fns = nii_fns(contains(nii_fns,{'P0252374'}));

% 
%-----------------------------------------
% Make paths to nii_xps, boostraps, and maps folders

method = 'dtd';

pmaps_paths = cell(0,0);
bs_paths = cell(0,0);
for ndata = 1:numel(nii_fns)
    rwi_path = fileparts(fileparts(nii_fns{ndata}));
    pmaps_paths{1+numel(pmaps_paths)} = fullfile(rwi_path,'pmaps');
    bs_paths{1+numel(bs_paths)} = fullfile(fileparts(fileparts(nii_fns{ndata})),method,'bootstraps');        
end

% Define bins
% disomin = [0 0 2.5 0 0]*1e-9; disomax = [2.5 2.5 5 5 5]*1e-9;
% sddeltamin = [.25 0 0 .25 0]; sddeltamax = [1 .25 1 1 .25];
disomin = [0 0 2.5]*1e-9; disomax = [2.5 2.5 5]*1e-9;
sddeltamin = [.25 0 0]; sddeltamax = [1 .25 1];
dratiomin = ones(size(disomin))*eps; dratiomax = ones(size(disomin))/eps;
r2min = .5*ones(size(disomin)); r2max = 30*ones(size(disomin));
r1min = .1*ones(size(disomin)); r1max = 2*ones(size(disomin));

% Define color limits for parameter maps
clim.s0 = 1*[0 1]; % Multiplied with s0 below
clim.s2000 = 1*[0 1]; % Multiplied with s0 below
clim.mdiso = 3.5e-9*[0 1]; %clim.mdiso = 1e-9*[0 1];
clim.msddelta = 1*[0 1];
clim.mr2 = max(r2max)*[0 1];
clim.mr1 = max(r1max)*[0 1]; clim.mr1 = .8*[0 1];
clim.vdiso = .3*3e-9^2*[0 1]; %clim.vdiso = .3*1e-9^2*[0 1];
clim.vsddelta = .15*[0 1];
clim.vr2 = .2*max(r2max)^2*[0 1];
clim.vr1 = .2*max(r1max)^2*[0 1]; clim.vr1 = .2*.8^2*[0 1];
clim.cvdisosddelta = .1*3e-9*1*[-1 1]; %clim.cvdisosddelta = .2*1e-9*1*[-1 1];
clim.cvdisor2 = .1*max(r2max)*3e-9*[-1 1];
clim.cvsddeltar2 = .1*max(r2max)*1*[-1 1];
clim.cvdisor1 = .1*max(r1max)*3e-9*[-1 1]; clim.cvdisor1 = .1*.8*3e-9*[-1 1];
clim.cvsddeltar1 = .1*max(r1max)*1*[-1 1]; clim.cvsddeltar1 = .1*.8*1*[-1 1];
clim.mask_threshold = eps; clim.mask_threshold = 0.02;

%------------------------------

% Prepare options
opt.(method).bin_disomin = disomin; opt.(method).bin_disomax = disomax;
opt.(method).bin_dratiomin = dratiomin; opt.(method).bin_dratiomax = dratiomax;
opt.(method).bin_sddeltamin = sddeltamin; opt.(method).bin_sddeltamax = sddeltamax;
if strcmp(method,'dtr2d')
    opt.(method).bin_r2min = r2min; opt.(method).bin_r2max = r2max;
elseif strcmp(method,'dtr1d')
    opt.(method).bin_r1min = r1min; opt.(method).bin_r1max = r1max;
end
opt = mdm_opt(opt);


Ndata = numel(bs_paths);
tic

% Loop over datasets
for ndata = 1:Ndata
    
    bs_dps = mdm_dps_collectbs(method, bs_paths{ndata}, opt);
        
    if ~all(cellfun('isempty',bs_dps))
        median_dps = mdm_dps_median(bs_dps);
        clear bs_dps
%%  
% clim.mask_threshold = 0.15;
%          mplot_technicolor_nii(method, median_dps, pmaps_paths{ndata}, clim, opt)
%           mplot_technicolor(method, median_dps, fullfile(fileparts(pmaps_paths{ndata}),'technicolor'), clim);
%          mplot_technicolor_slices(method, median_dps, fullfile(fileparts(pmaps_paths{ndata}),'slices'), clim, opt)
         mplot_technicolor_slicecollage(method, median_dps, fullfile(fileparts(pmaps_paths{ndata}),'slicecollage'), clim, opt)
%           mplot_globalstats(method, median_dps, fullfile(fileparts(pmaps_paths{ndata}),'globalstats'), clim)
    end
   
end
toc


