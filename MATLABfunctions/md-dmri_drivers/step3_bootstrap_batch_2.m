%Run bootstrap analysis
clear all

wd = pwd;

method = 'dtd';

%Define paths to bootstrap folders
datasets_paths{1} = '/Users/daniel/Dropbox/NMRdata/NIH_Philips/mdd';
% 
Ndata = numel(datasets_paths);
in_paths = cell(0,0);
bs_paths = cell(0,0);
out_paths = cell(0,0);
for ndata = 1:Ndata
    datasets_path = datasets_paths{ndata};
    expnams = mdm_bruker_dir2expnams(datasets_path);
    Nexp = numel(expnams);
    for nexp = 1:Nexp
        expnam = expnams{nexp};
        in_paths{1+numel(in_paths)} = fullfile(datasets_path,expnam,'nii_xps');
        bs_paths{1+numel(bs_paths)} = fullfile(datasets_path,expnam,method,'bootstraps');
        out_paths{1+numel(out_paths)} = fullfile(datasets_path,expnam,method,'maps');
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

opt = dtd_opt(opt);
opt.dtd.ind_start = 1;
opt.dtd.dmin = 5e-11;
opt.dtd.dmax = 5e-9;
opt.dtd.n_out = 10;


for ndata = 1:Ndata
    in_path = in_paths{ndata};
    bs_path = bs_paths{ndata};

    % Connect to data
    clear s
    s.nii_fn = fullfile(in_path, [nii_name '.nii.gz']);
    s.mask_fn = fullfile(in_path, [nii_name '_mask.nii.gz']);
    s.xps = mdm_xps_load(fullfile(in_path, [nii_name '_xps.mat']));
    
    opt.mask.b0_ind = find(abs(s.xps.b-min(s.xps.b)) < 1e7);

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

