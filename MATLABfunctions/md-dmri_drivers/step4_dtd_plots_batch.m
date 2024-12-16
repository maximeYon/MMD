%Plot parameters maps and global stats
clear all

% Define paths to dataset folders
datasets_paths = cell(0);
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/USC/Data3';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients_20200130';
datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/United/data4dtdpaper';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd19';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd18';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd17';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/USC/Data1';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd16';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd15';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd14';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd13';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd12';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd11';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd10';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd09';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd08';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd07';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd06';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd05';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd04';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd03';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd02';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191029night_v1';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191029pm_v1';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/GE/2018-11-22_Premier_corrected_waveforms';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/NIH_Philips';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/MR750_20190521/RWI_volunteer';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191022';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191023';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191023_pm';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Caeyenberghs/FWF_DTC01_Caeyenberghs';
% 
%-----------------------------------------
% Make paths to nii_xps, boostraps, and maps folders

method = 'dtd';

paths.nii_xps = cell(0,0);
paths.bootstraps = cell(0,0);
paths.maps = cell(0,0);
for ndata = 1:numel(datasets_paths)
    datasets_path = datasets_paths{ndata};
    mdd_path = fullfile(datasets_path,'mdd');
    nifti_names = mdm_bruker_dir2expnams(mdd_path);
    
    for nnifti = 1:numel(nifti_names)
        nifti_name = nifti_names{nnifti};
        nifti_path = fullfile(mdd_path, nifti_name); 
        paths.nii_xps{1+numel(paths.nii_xps)} = fullfile(nifti_path,'nii_xps');
        paths.bootstraps{1+numel(paths.bootstraps)} = fullfile(nifti_path,method,'bootstraps');
        paths.maps{1+numel(paths.maps)} = fullfile(nifti_path,method,'maps');
    end
end

% Define bins
disomin = [0 0 2.5]*1e-9; disomax = [2.5 2.5 5]*1e-9;
dratiomin = [1 1 1]*eps; dratiomax = [1 1 1]/eps;
sddeltamin = [.25 0 0]; sddeltamax = [1 .25 1];
% disomin = [0 0 2]*1e-9; disomax = [2 2 5]*1e-9;
% dratiomin = [1 1 1]*eps; dratiomax = [1 1 1]/eps;
% sddeltamin = [.25 0 0]; sddeltamax = [1 .25 1];
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

% Define color limits for parameter maps
clim.s0 = .5*[0 1]; % Multiplied with s0 below
clim.s2000 = .5*[0 1]; % Multiplied with s0 below
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

% Prepare options
opt = mdm_opt();
opt.(method).bin_disomin = disomin; opt.(method).bin_disomax = disomax;
opt.(method).bin_dratiomin = dratiomin; opt.(method).bin_dratiomax = dratiomax;
opt.(method).bin_sddeltamin = sddeltamin; opt.(method).bin_sddeltamax = sddeltamax;
if strcmp(method,'dtr2d')
    opt.(method).bin_r2min = r2min; opt.(method).bin_r2max = r2max;
elseif strcmp(method,'dtr1d')
    opt.(method).bin_r1min = r1min; opt.(method).bin_r1max = r1max;
end

Ndata = numel(paths.nii_xps);
tic

% Loop over datasets
for ndata = 1:Ndata

    bs_dps = mdm_dps_collectbs(method, paths.bootstraps{ndata}, opt);
        
    if ~all(cellfun('isempty',bs_dps))
        median_dps = mdm_dps_median(bs_dps);
        clear bs_dps
       
         mplot_technicolor_nii(method, median_dps, fullfile(paths.maps{ndata},'nii'), clim, opt)
         mplot_technicolor(method, median_dps, fullfile(paths.maps{ndata},'technicolor'), clim, opt)
 %         mplot_technicolor_slices(method, median_dps, fullfile(paths.maps{ndata},'slices'), clim, opt)
%          mplot_technicolor_slicecollage(method, median_dps, fullfile(paths.maps{ndata},'slicecollage'), clim, opt)
%         mplot_globalstats(method, median_dps, fullfile(paths.maps{ndata},'globalstats'), clim, opt)
    end
   
end
toc

