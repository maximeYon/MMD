%Run bootstrap analysis
clear all

% Define paths to dataset folders
% data_path =  '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419';
data_path =  '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
data_dir = dir(fullfile(data_path,'*'));

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'2020*'));

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'20*'));

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_Mahajan_1/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'20*'));

% data_path = '/Users/daniel/Dropbox/NMRdata/Philips_QC/acquisition_date_time';
% % data_path = '/Users/daniel/Dropbox/NMRdata/GE_QC/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'*'));

nii_fns = cell(size(data_dir,1),1);
for ndata = 1:size(data_dir,1)
    if ~contains(data_dir(ndata).name,'.')
        nii_fns{ndata} = fullfile(data_path,data_dir(ndata).name,'rwi','nii_xps','data_mc.nii.gz');
    end
end
nii_fns = nii_fns(~cellfun(@isempty,nii_fns));
% nii_fns = nii_fns(end:-1:1);
% nii_fns = nii_fns(contains(nii_fns,'p0235177'));

%-----------------------------------------
% Select fit methods
methods = {'dti_lls'; 'dtd_covariance'; 'dti_euler'; 'dtd_gamma'};
methods = {'dtd_gamma'};

% Define paths to parameter map folders
pmaps_paths = cell(0,0);
for ndata = 1:numel(nii_fns)
    rwi_path = fileparts(fileparts(nii_fns{ndata}));
    pmaps_paths{1+numel(pmaps_paths)} = fullfile(rwi_path,'pmaps');
end

% Prepare options
opt = mdm_opt();
opt = mio_opt(opt);
opt.verbose       = 1;
opt.do_overwrite  = 1;
opt.do_mask = 0; 
opt.mask.do_overwrite = 0;
opt.mask.threshold = .1;
opt.do_data2fit = 0; 
opt.do_fit2param = 1;
opt.do_param2nii = 1;
opt.mio.do_parfor = 1;

opt = dti_lls_opt(opt);

opt = dtd_covariance_opt(opt);
opt.dtd_covariance.do_clamping = 1;
opt.dtd_covariance.fig_maps = {'s0','s1000','s2000','MD','FA','nmsdaniso','nvdiso','MKi','MKa','MKt'};

opt = dti_euler_opt(opt);
opt.dti_euler.weight_sthresh = 0.1;

opt = dtd_gamma_opt(opt);
opt.dtd_gamma.weight_sthresh = 0.1;
opt.dtd_gamma.do_multiple_s0 = 0;
opt.dtd_gamma.do_pa = 1;
opt.dtd_gamma.do_clamping = 1;
opt.dtd_gamma.fig_maps = {'s0','s1000','s2000','MD','nmsdaniso','nvdiso','MKi','MKa','MKt'};

global_opt = opt;

% Loop over datasets
for ndata = 1:numel(nii_fns)
    % Connect to data
    clear s
    s.nii_fn = nii_fns{ndata};
    s.mask_fn = mdm_fn_nii2mask(s.nii_fn, opt);
    s.xps = mdm_xps_load(mdm_fn_nii2xps(s.nii_fn));
        
    opt.mask.b0_ind = find(s.xps.b < .11e9);    
        
    %Loop over methods
    for nmethod = 1:numel(methods)
        method = methods{nmethod};

        method_path = fullfile(fileparts(pmaps_paths{ndata}),method);
        msf_mkdir(method_path);

        clear paths
        paths.nii_path   = pmaps_paths{ndata};
        paths.mfs_fn   = fullfile(method_path, 'mfs.mat');
        paths.dps_fn   = fullfile(method_path, 'dps.mat');
        paths = mdm_paths(paths);
        
        % Save opt for book keeping
        opt_fn = fullfile(method_path,'opt.mat');
        if opt.do_overwrite == 1 || exist(opt_fn,'file')~=2      
            opt = global_opt;
            save(opt_fn,'opt')
        end

        % Run analysis
        tic;
        pipename = [method '_pipe'];
        nii_fn = feval(pipename, s, paths, opt);
        toc;
    end
end



