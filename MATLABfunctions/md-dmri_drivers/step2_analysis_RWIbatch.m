%Run bootstrap analysis
clear all

% Define paths to dataset folders
data_paths = cell(0);
data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0240306/rwi/res_mc';
data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0275752/rwi/res_mc';
data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0333068/rwi/res_mc';

%-----------------------------------------
% Make paths to nii_xps, boostraps, and maps folders

method = 'dti_euler';
nii_name = 'data_mc';

paths.nii_xps = cell(0,0);
paths.maps = cell(0,0);
for ndata = 1:numel(data_paths)
    data_path = data_paths{ndata};
    paths.nii_xps{1+numel(paths.nii_xps)} = data_path;
    paths.maps{1+numel(paths.maps)} = fullfile(fileparts((data_path)),['res_' method]);
end

% Prepare options
opt = mdm_opt();
opt = mio_opt(opt);
opt.verbose       = 1;
opt.do_overwrite  = 1;
opt.do_mask = 1; 
opt.mask.do_overwrite = 1;
opt.mask.threshold = .1;
opt.do_data2fit = 1; 
opt.do_fit2param = 1;
opt.mio.do_parfor = 1;

opt = dti_euler_opt(opt);
opt.dti_euler.weight_sthresh = 0.1;

% Loop over datasets
for ndata = 1:numel(paths.nii_xps)

    % Connect to data
    clear s
    s.nii_fn = fullfile(paths.nii_xps{ndata}, [nii_name '.nii.gz']);
    s.mask_fn = fullfile(paths.nii_xps{ndata}, [nii_name '_mask.nii.gz']);
    s.xps = mdm_xps_load(fullfile(paths.nii_xps{ndata}, [nii_name '_xps.mat']));
        
    opt.mask.b0_ind = find(abs(s.xps.b-min(s.xps.b)) < .11e9);    
        
    % Run analysis
    tic;

    o     = fullfile(paths.maps{ndata},[method '_res1']);
    msf_mkdir(o);
    opt_fn = fullfile(o,'opt.mat');
    if exist(opt_fn,'file')==2            
        tmp = load(opt_fn,'opt'); opt = tmp.opt;
    else
        save(opt_fn,'opt')
    end

    paths.nii_path   = o;
    paths.mfs_fn   = fullfile(o, [method '_mfs.mat']);
    paths.dps_fn   = fullfile(o, [method '_dps.mat']);
    paths = mdm_paths(paths);

    nii_fn = dti_euler_pipe(s, paths, opt);


    toc;
end

