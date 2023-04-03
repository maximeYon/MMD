%Run bootstrap analysis
clear all

% Define full paths to the nifti files to be analyzed 
data_path =  '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
data_dir = dir(fullfile(data_path,'P025237*'));

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'2020*'));

nii_fns = cell(size(data_dir,1),1);
for ndata = 1:size(data_dir,1)
    nii_fns{ndata} = fullfile(data_path,data_dir(ndata).name,'rwi','nii_xps','data_mc.nii.gz');
end
% nii_fns = nii_fns(end:-1:1);



% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191002090825/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191002184711/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191010133317/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191010133850/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191021111255/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191024170921/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191024171722/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191025191039/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191107094708/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191108140810/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191111084057/rwi/nii_xps';
% data_paths{1+numel(data_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20200121173005/rwi/nii_xps';

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


% Prepare options
opt = mdm_opt();
opt.do_mask      = 0;
opt.mask.threshold = 0.02;
opt.mio.no_parfor = 1;
opt.do_data2fit  = 1;
opt.do_bootstrap = 1;
opt.do_fit2param = 0;
opt.do_param2nii = 0;

opt.do_xps2pdf    = 0;
opt.do_nii2pdf    = 0;
opt.do_m2pdf      = 0;

opt.verbose       = 1;
opt.do_overwrite  = 1;

if strcmp(method,'dtr2d')
    opt = dtr2d_opt(opt);
    opt.dtr2d.ind_start = 1;
    opt.dtr2d.dmin = 5e-11;
    opt.dtr2d.dmax = 5e-9;
    opt.dtr2d.r2min = 1;
    opt.dtr2d.r2max = 30;
    opt.dtr2d.n_out = 10;
elseif strcmp(method,'dtr1d')
    opt = dtr1d_opt(opt);
    opt.dtr1d.ind_start = 1;
    opt.dtr1d.dmin = 5e-11;
    opt.dtr1d.dmax = 5e-9;
    opt.dtr1d.r1min = .1;
    opt.dtr1d.r1max = 1;
    opt.dtr1d.n_out = 10;
elseif strcmp(method,'dtd')
    opt = dtd_opt(opt);
    opt.dtd.ind_start = 1;
    opt.dtd.dmin = 5e-11;
    opt.dtd.dmax = 5e-9;
    opt.dtd.n_out = 15;
end

% Loop over datasets
for ndata = 1:numel(nii_fns)
    % Connect to data
    clear s
    s.nii_fn = nii_fns{ndata};
    s.mask_fn = mdm_fn_nii2mask(s.nii_fn, opt);
    s.xps = mdm_xps_load(mdm_fn_nii2xps(s.nii_fn));
            
    bs_path = bs_paths{ndata};

    msf_mkdir(bs_path);
    opt_fn = fullfile(bs_path,'opt.mat');
    if exist(opt_fn,'file')==2            
        tmp = load(opt_fn,'opt'); opt = tmp.opt;
    else
        save(opt_fn,'opt')
    end
    
    % Run analysis
    tic;

    parfor nbs = 1:96
        o     = fullfile(bs_path,num2str(nbs));
        msf_mkdir(o);

        paths.nii_path   = fullfile(o, 'nii_res');
        paths.mfs_fn   = fullfile(o, 'mfs.mat');
        paths.ind_fn   = fullfile(o, 'ind.mat');
        paths.dps_fn   = fullfile(o, 'dps.mat');
        paths = mdm_paths(paths);
                
        if exist(paths.mfs_fn,'file')~=2            
            pipename = [method '_pipe'];
            nii_fn = feval(pipename, s, paths, opt);
        end

    end

    toc;
end

