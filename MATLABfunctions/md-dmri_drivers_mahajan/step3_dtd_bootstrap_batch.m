%Run bootstrap analysis
clear all

% Define paths to dataset folders
datasets_paths = cell(0);
datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd16';
datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd15';
datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd14';
datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd13';
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

nii_name = 'data';

% Prepare options
opt = mdm_opt();
opt.do_mask      = 0;
opt.mask.threshold = 0.1;
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
    opt.dtr2d.dmin = 5e-12;
    opt.dtr2d.dmax = 5e-9;
    opt.dtr2d.r2min = 1;
    opt.dtr2d.r2max = 80;
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
    opt.dtd.n_out = 10;
elseif strcmp(method,'dti_euler')
    opt = dti_euler_opt(opt);
    opt.dti_euler_opt.ind_start = 1;
end

% Loop over datasets
for ndata = 1:numel(paths.nii_xps)

    % Connect to data
    clear s
    s.nii_fn = fullfile(paths.nii_xps{ndata}, [nii_name '.nii.gz']);
    s.mask_fn = fullfile(paths.nii_xps{ndata}, [nii_name '_mask.nii.gz']);
    s.xps = mdm_xps_load(fullfile(paths.nii_xps{ndata}, [nii_name '_xps.mat']));
        
    msf_mkdir(paths.bootstraps{ndata});
    opt_fn = fullfile(paths.bootstraps{ndata},'opt.mat');
    if exist(opt_fn,'file')==2            
        tmp = load(opt_fn,'opt'); opt = tmp.opt;
    else
        save(opt_fn,'opt')
    end
    
    % Run analysis
    tic;

    path_bootstraps = paths.bootstraps{ndata};
    parfor nbs = 1:96
        o     = fullfile(path_bootstraps,num2str(nbs));
        msf_mkdir(o);

        paths.nii_path   = fullfile(o, 'nii_res');
        paths.mfs_fn   = fullfile(o, 'mfs.mat');
        paths.ind_fn   = fullfile(o, 'ind.mat');
        paths.dps_fn   = fullfile(o, 'dps.mat');
        paths = mdm_paths(paths);
        
        if exist(paths.mfs_fn,'file')~=2            
            if strcmp(method,'dtr2d')
                nii_fn = dtr2d_pipe(s, paths, opt);
            elseif strcmp(method,'dtr1d')
                nii_fn = dtr1d_pipe(s, paths, opt);
            elseif strcmp(method,'dtd')
                nii_fn = dtd_pipe(s, paths, opt);
            elseif strcmp(method,'dti_euler')
                nii_fn = dti_euler_pipe(s, paths, opt);
            end
        end

    end

    toc;
end

