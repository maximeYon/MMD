%Run bootstrap analysis
clear all

wd = pwd;

%method = 'dti_euler';
method = 'dtd';
%method = 'dtr1d';

%Define paths to bootstrap folders
% datasets_path = '/Users/daniel/Dropbox/NMRdata/United/data4dtdpaper';
% datasets_path = '/Users/daniel/Dropbox/NMRdata/Caeyenberghs/FWF_DTC01_Caeyenberghs/mdd';
% datasets_path = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd';
% expnams = mdm_bruker_dir2expnams(datasets_path);
% %expnams = expnams(4);
% 
% Ndata = numel(expnams);
% bs_paths = cell(Ndata,1);
% out_paths = cell(Ndata,1);
% for ndata = 1:Ndata
%     expnam = expnams{ndata};
%     in_paths{ndata} = fullfile(datasets_path,expnam,'nii_xps');
%     bs_paths{ndata} = fullfile(datasets_path,expnam,method,'bootstraps');
% end


% datasets_path = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd';
% expnams = mdm_bruker_dir2expnams(datasets_path);
% expnams = expnams(3);
% 
% Ndata = numel(expnams);
% in_paths = cell(Ndata,1);
% bs_paths = cell(Ndata,1);
% out_paths = cell(Ndata,1);
% for ndata = 1:Ndata
%     expnam = expnams{ndata};
%     in_paths{ndata} = fullfile(datasets_path,expnam,'mc','nii_xps');
%     bs_paths{ndata} = fullfile(datasets_path,expnam,'mc',method,'bootstraps');
% end

% 
%Define paths to bootstrap folders
datasets_paths{1} = '/Users/daniel/Dropbox/NMRdata/GE/2018-11-22_Premier_corrected_waveforms/mdd';
%datasets_paths{1} = '/Users/daniel/Dropbox/NMRdata/NIH_Philips/mdd';
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

% 
% clear in_paths bs_paths
% method = 'dtr2d';
% in_paths{1} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_15_20190725_161749/mdd/nii_xps_merge_21to30';
% in_paths{1+numel(in_paths)} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_17_20190726_115920/mdd/nii_xps_merge_21to30';
% in_paths{1+numel(in_paths)} = '/Users/daniel/Dropbox/NMRdata/MaglabTallahassee/July2019/Topgaard_D_1_14_20190725_122035/mdd/nii_xps_merge';
% Ndata = numel(in_paths);
% for ndata = 1:Ndata
%     bs_paths{ndata} = fullfile(in_paths{ndata},method,'bootstraps');
% end


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

Ndata = numel(in_paths);
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

