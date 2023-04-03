%Run bootstrap analysis
clearvars;
parpool();

home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

% Define full paths to the nifti files to be analyzed 
%data_path = '20211019_110649_test_sample_Hong_test_sample_Hong_1_1\127';
% data_path = '20220208_HongPhantom3\9';
% data_path = '20220204_tumor3\60'; 
% data_path = '20220202_RatBrain\22'; 
% data_path = 'Original_Rat_Brain\14'; 
% data_path = '20220407_Mouse_heart_PBS\65';
%data_path = '20220411_Mouse_Brain_1\32'; 
% data_path = '20220420_Mouse_brain2\15'; 
% data_path = '20220524_mouse_brain_fin2\70';
%data_path = 'test\321'; 
% data_path = '20220620_130815_qMAS_testing\132'; 
% data_path = '20220722_mouse_brain_fin4\8'; 
% data_path = 'spectro\37';
% data_path = 'omar_hippocampus_sample\10'; 
% data_path = '20221206_kenneth_data\30'; 
% data_path = '20221207_fruit\48';
%data_path = '20221211_mouse_heart\22';
% data_path = '\20221209_Hong_phantom_spectro\5';
% data_path = 'clay_1\13';
data_path = 'human_sample_fresh1\34';


data_path = split(data_path,'\'); data_path = join(data_path,filesep,1); data_path = data_path{1};
% nii_fns{1} = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep 'pdata_mdd' filesep 'nii_xps' filesep 'data.nii'];
nii_fns{1} = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep 'pdata_mdd' filesep 'nii_xps' filesep 'data.nii.gz'];

%-----------------------------------------
% Make paths to nii_xps, boostraps, and maps folders

% method = 'dtd';
% method = 'dtod';
method = 'dtor1r2d';
% method = 'dtr1r2d';

pmaps_paths = cell(0,0);
bs_paths = cell(0,0);
for ndata = 1:numel(nii_fns)
    pdata_path = fileparts(fileparts(nii_fns{ndata}));
    pmaps_paths{1+numel(pmaps_paths)} = fullfile(pdata_path,'pmaps');
    bs_paths{1+numel(bs_paths)} = fullfile(pdata_path,method,'bootstraps');        
end

bsno = 1:100; % Bootstrap numbers

% Prepare options
opt = mdm_opt();
opt = feval([method '_opt'],opt);
opt.do_mask      = 1;
opt.mask.threshold = 0.1;
opt.mio.do_parfor = 1;
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
    opt.dtr2d.ind_start = 1;
    opt.dtr2d.dmin = 5e-11;
    opt.dtr2d.dmax = 5e-9;
    opt.dtr2d.r2min = 1;
    opt.dtr2d.r2max = 30;
    opt.dtr2d.n_out = 10;
elseif strcmp(method,'dtr1d')
    opt.dtr1d.ind_start = 1;
    opt.dtr1d.dmin = 5e-11;
    opt.dtr1d.dmax = 5e-9;
    opt.dtr1d.r1min = .1;
    opt.dtr1d.r1max = 1;
    opt.dtr1d.n_out = 10;
elseif strcmp(method,'dtr1r2d')
    opt.(method).ind_start = 1;
    opt.(method).dmin = 5e-12;
    opt.(method).dmax = 5e-9;
    opt.(method).r1min = .1;
    opt.(method).r1max = 4;
    opt.(method).r2min = 4;
    opt.(method).r2max = 150;
    opt.(method).n_out = 10;
elseif strcmp(method,'dtor1r2d')
    opt.(method).ind_start = 1;
    opt.(method).dmin = 5e-12;
    opt.(method).dmax = 5e-9;
    opt.(method).rmin = .1;
    opt.(method).rmax = 1e5;
    opt.(method).r1min = .1;
    opt.(method).r1max = 4;
    opt.(method).r2min = 4; 
    opt.(method).r2max = 150; %150
    opt.(method).n_out = 10;
elseif strcmp(method,'dtod')
    opt.(method).ind_start = 1;
    opt.(method).dmin = 5e-12;% min diffusion limits
    opt.(method).dmax = 5e-9;% max diffusion limits
    opt.(method).rmin = .1;% min restriction limits
    opt.(method).rmax = 1e5;% max restriction limits
    opt.(method).n_out = 10;% Number of output components
elseif strcmp(method,'dtd')
    opt.dtd.ind_start = 1;
    opt.dtd.dmin = 5e-12;
    opt.dtd.dmax = 5e-9;
    opt.dtd.n_out = 10;
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

    for nbs = bsno
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

