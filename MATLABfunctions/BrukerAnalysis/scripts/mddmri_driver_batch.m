% Driver for processing data acquired with Bruker pulse sequence DT_axderare2d.% Topgaard, Phys. Chem. Chem. Phys. (2016).% http://dx.doi.org/10.1039/c5cp07251d.% Assumes axisymmetric b-tensors and powder averaging.% % Select models% 1) Conventional diffusion tensors;%    fractional anisotropy;%    Westin's shape indices.%% 2) Isotropic and anisotropic variance of the diffusion tensor distribution;%    orientational order parameters;%    microscopic diffusion anisotropy.%    See Lasic et al, Front. Phys. 2, 11 (2014). %    http://dx.doi.org/10.3389/fphy.2014.00011.%% 3) Shape of the microscopic diffusion tensor (prolate, sphere, oblate).%    See, Eriksson et al., J. Chem. Phys. 142, 104201 (2015).%    http://dx.doi.org/10.1063/1.4913502.%% 4) Saupe order tensors.%    See Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016).%    http://dx.doi.org/10.1039/c5cp07251d.%% 5) Size-shape diffusion tensor distributions.%    See de Almeida Martins and Topgaard, Phys. Rev. Lett. 116, 087601 (2016).%    http://dx.doi.org/10.1103/PhysRevLett.116.087601.%   % 6) Size-shape-orientation diffusion tensor distributions.%    See Topgaard. J. Magn. Reson. 275, 98 (2017).%    http://dx.doi.org/10.1016/j.jmr.2016.12.007clear allmodels         = {'dti_euler', 'dtd_gamma', 'dtd_pake', 'dtd_saupe', 'dtd_pa', 'dtd'};c_model        = 6;% Define path for data and output foldersdatadir = '/Users/daniel/Dropbox/NMRdata/AVII500';nnam = 0;%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq17'; datacell{nnam}.no = [3:188]; datacell{nnam}.no = [155:188];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq15'; datacell{nnam}.no = [5:204];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq16'; datacell{nnam}.no = [5:304];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq14'; datacell{nnam}.no = [5:217]; datacell{nnam}.no = [3:4];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq13'; datacell{nnam}.no = [7:75]; datacell{nnam}.no = 75;%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq12'; datacell{nnam}.no = [5:190];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq11'; datacell{nnam}.no = [5:204];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq10'; datacell{nnam}.no = [3:202];%nnam = nnam+1; datacell{nnam}.nam = 'C14E5_Eq'; datacell{nnam}.no = [3 8:107]; datacell{nnam}.no = [108:206];%nnam = nnam+1; datacell{nnam}.nam = 'DTD2'; datacell{nnam}.no = [3:14 16:20 22 24:46]; datacell{nnam}.no = [47:131];%nnam = nnam+1; datacell{nnam}.nam = 'DT_axderare2d_test'; datacell{nnam}.no = [4 7 10 13]; %datacell{nnam}.no = [13];%nnam = nnam+1; datacell{nnam}.nam = 'DTD'; datacell{nnam}.no = [11:121]; datacell{nnam}.no = [21];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_temp12'; datacell{nnam}.no = [2:45];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_temp13'; datacell{nnam}.no = [1:84];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq8'; datacell{nnam}.no = [98:100];%nnam = nnam+1; datacell{nnam}.nam =  'AOToct_temp14'; datacell{nnam}.no = [17:26];%nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq9'; datacell{nnam}.no = [38:150];%nnam = nnam+1; datacell{nnam}.nam = 'MIC5_2H_setup'; datacell{nnam}.no = []; datacell{nnam}.no = [66:67];nnam = nnam+1; datacell{nnam}.nam = 'AOToct_Eq18'; datacell{nnam}.no = [3 6:213];% Image recon parametersrps.smooth          = 300e-6; %Gaussian smoothing in [m]rps.npix.read       = 16; % Number of pixels in read dimensionrps.npix.phase      = rps.npix.read; % Number of pixels in phase dimensionrps.shift_read      = -.25e-3; % Image shift in read dimension [m]% Prepare optionsopt = mdm_opt();opt.do_recon     = 1;opt.do_mask      = 1;opt.do_data2fit  = 1;opt.do_fit2param = 1;opt.do_param2nii = 1;opt.do_xps2pdf    = 0;opt.do_nii2pdf    = 0;opt.do_dtdpdf     = 1;opt.verbose       = 1;opt.do_overwrite  = 1;opt.dti_euler.fig_maps  = {'s0','iso','fa'};opt.dtd_gamma.fig_maps      = {'s0','iso','ciso','cmu'};opt.dtd_pake.fig_maps        = {'s0','iso','delta'};opt.dtd.fig_maps        = {};opt.dtd.fig_cmaps        = {};opt.mask.thresh         = 0.05;for nnam = 1:numel(datacell)    expnam = datacell{nnam}.nam;    expno = datacell{nnam}.no;    for nno = 1:numel(expno)        data_path = fullfile(datadir,expnam,num2str(expno(nno)));        out_path  = fullfile(data_path, 'NII_XPS');         msf_mkdir(out_path);        if (opt.do_recon)           mdm_bruker_dt_rare2d_recon(data_path, out_path, rps);        end                i     = fullfile(data_path, 'NII_XPS');         o     = fullfile(data_path, 'NII_RES');        msf_mkdir(o);        % Connect to data        s.nii_fn = fullfile(i, 'data_sub.nii.gz');        s.xps = mdm_xps_load(fullfile(i, 'data_sub_xps.mat'));        if (opt.do_xps2pdf)           mdm_xps2pdf(i,opt);        end        % Run analysis        for n_model = 1:numel(c_model)            tic;            % OUTPUT: define paths for data, fit parameters, and maps            paths.nii_path = fullfile(o, 'maps');            paths.mfs_fn   = fullfile(o, models{c_model(n_model)}, 'mfs.mat');            paths.dps_fn   = fullfile(o, models{c_model(n_model)}, 'dps.mat');            switch (models{c_model(n_model)})                case 'dti_euler'                    nii_fn = dti_euler_pipe(s, paths, opt);                case 'dtd_gamma'                    nii_fn = dtd_gamma_pipe(s, paths, opt);                case 'dtd_pake'                    nii_fn = dtd_pake_pipe(s, paths, opt);                case 'dtd_saupe'                    nii_fn = dtd_saupe_pipe(s, paths, opt);                case 'dtd_pa'                    nii_fn = dtd_pa_pipe(s, paths, opt);                case 'dtd'                    nii_fn = dtd_pipe(s, paths, opt);            end            toc;        end    endend