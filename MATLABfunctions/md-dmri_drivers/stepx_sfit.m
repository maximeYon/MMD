%Run bootstrap analysis
clear all

wd = pwd;

%method = 'dti_euler';
method = 'dtd';
%method = 'dtr1d';

%Define paths to bootstrap folders
% datasets_path = '/Users/daniel/Dropbox/NMRdata/United/data4dtdpaper';
% expnams = mdm_bruker_dir2expnams(datasets_path);
% expnams = expnams(1);
% 
% Ndata = numel(expnams);
% bs_paths = cell(Ndata,1);
% out_paths = cell(Ndata,1);
% for ndata = 1:Ndata
%     expnam = expnams{ndata};
%     in_paths{ndata} = fullfile(datasets_path,expnam,'nii_xps');
%     bs_paths{ndata} = fullfile(datasets_path,expnam,method,'bootstraps');
% end

datasets_path = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191022/mdd';
datasets_path = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191023/mdd';
datasets_path = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191023_pm/mdd';
datasets_path = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191029pm_v1/mdd';
datasets_path = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191029night_v1/mdd';
datasets_path = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191030am_p1/mdd';
datasets_path = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd';
datasets_path = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191109_p04/mdd';

expnams = mdm_bruker_dir2expnams(datasets_path);
%expnams = expnams([3])
%expnams = expnams(fliplr(1:numel(expnams)))
%return
Ndata = numel(expnams);
bs_paths = cell(Ndata,1);
out_paths = cell(Ndata,1);
for ndata = 1:Ndata
    expnam = expnams{ndata};
    in_paths{ndata} = fullfile(datasets_path,expnam,'nii_xps');
    bs_paths{ndata} = fullfile(datasets_path,expnam,method,'bootstraps');
end

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

for ndata = 1:Ndata
    in_path = in_paths{ndata};
    bs_path = bs_paths{ndata};

    % Connect to data
    clear s
    s.nii_fn = fullfile(in_path, [nii_name '.nii.gz']);
    s.mask_fn = fullfile(in_path, [nii_name '_mask.nii.gz']);
    s.xps = mdm_xps_load(fullfile(in_path, [nii_name '_xps.mat']));
    
    [I,h]   = mdm_nii_read(s.nii_fn);
    M = mdm_mask_load(s, opt);
    xps = s.xps;

    for nbs = 1:96
        o     = fullfile(bs_path,num2str(nbs));
        paths.mfs_fn   = fullfile(o, 'mfs.mat');
        if exist(paths.mfs_fn,'file')==2

            mfs = mdm_mfs_load(paths.mfs_fn);

            m = mfs.m;        
            sz = size(m);

            I_fit = zeros(sz(1),sz(2),sz(3),xps.n);
            for nk = 1:sz(3)
                for nj = 1:sz(2)
                    for ni = 1:sz(1)
                        if M(ni,nj,nk)
                            %s_fit = dtr2d_1d_fit2data(squeeze(m(ni,nj,nk,:))', xps);
                            s_fit = dtd_1d_fit2data(squeeze(m(ni,nj,nk,:))', xps);
                            I_fit(ni,nj,nk,:) = s_fit;
                            %figure(1), clf, plot(1:xps.n,squeeze(I_fit(ni,nj,nk,:)),'.'), pause(.1)
                        end
                    end
                end
            end

            nii_fit_fn   = fullfile(fileparts(s.nii_fn), ['data_fit_bs' num2str(nbs) '.nii.gz']);
            mdm_nii_write(I_fit, nii_fit_fn, h, 0);  
            xps_fit_fn   = fullfile(fileparts(s.nii_fn), ['data_fit_bs' num2str(nbs) '_xps.mat']);
            save(xps_fit_fn,'xps');

        end
    end
end
