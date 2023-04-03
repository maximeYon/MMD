%Plot parameters maps and global stats
clear all

% Define paths to dataset folders
data_paths = cell(0);
data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0240306/rwi/res_mc';
data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0275752/rwi/res_mc';
data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0333068/rwi/res_mc';

% 
%-----------------------------------------
% Make paths to nii_xps, boostraps, and maps folders

method = 'dtd';
nii_name = 'data_mc';

paths.nii_xps = cell(0,0);
paths.maps = cell(0,0);
paths.bootstraps = cell(0,0);
for ndata = 1:numel(data_paths)
    data_path = data_paths{ndata};
    paths.nii_xps{1+numel(paths.nii_xps)} = data_path;
    paths.maps{1+numel(paths.maps)} = fullfile(fileparts((data_path)),['res_' method],[method '_res1']);
    paths.bootstraps{1+numel(paths.bootstraps)} = fullfile(paths.maps{numel(paths.maps)},'bootstraps');        
end


% paths.nii_xps = cell(0,0);
% paths.nii_name = cell(0,0);
% paths.bootstraps = cell(0,0);
% paths.maps = cell(0,0);
% for ndata = 1:numel(datasets_paths)
%     datasets_path = datasets_paths{ndata};
%     %mdd_path = fullfile(datasets_path,'mdd');
%     mdd_path = datasets_path;
%     nifti_names = mdm_bruker_dir2expnams(mdd_path);    
%     for nnifti = 1:numel(nifti_names)
%         nifti_name = nifti_names{nnifti};
%         nifti_path = fullfile(mdd_path, nifti_name); 
%         nii_name = [nifti_name '_sub_mc'];
%         paths.nii_xps{1+numel(paths.nii_xps)} = fullfile(nifti_path);
%         paths.bootstraps{1+numel(paths.bootstraps)} = fullfile(fileparts(fileparts((nifti_path))),['res_' method],nifti_name,'bootstraps');
%         paths.maps{1+numel(paths.maps)} = fullfile(fileparts(fileparts((nifti_path))),['res_' method],nifti_name,'maps');
%         paths.nii_name{1+numel(paths.nii_name)} = nii_name;
%     end
% end

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
%%        
         mplot_technicolor_nii(method, median_dps, paths.maps{ndata}, clim, opt)
         mplot_technicolor(method, median_dps, fullfile(paths.maps{ndata},'technicolor'), clim)
% %         mplot_technicolor_slices(method, median_dps, fullfile(paths.maps{ndata},'slices'), clim)
% %          mplot_technicolor_slicecollage(method, median_dps, fullfile(paths.maps{ndata},'slicecollage'), clim)
        mplot_globalstats(method, median_dps, fullfile(paths.maps{ndata},'globalstats'), clim)
    end
   
end
toc


