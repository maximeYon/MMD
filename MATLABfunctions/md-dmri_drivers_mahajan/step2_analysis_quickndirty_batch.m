%    Quick n' dirty initial analysis of MDD data
%
%    Conventional per-voxel diffusion tensors
%
%    Mean-square anisotropy ("shape") and variance of isotropic
%    diffusivities ("size" heterogeneity) 
%
%    Parameters defined in Lasic et al, Front. Phys. 2, 11 (2014).
%    http://dx.doi.org/10.3389/fphy.2014.00011
%   
%    Quick covariance processing from Westin et al, Neuroimage 135, 345-362 (2016).
%    http://dx.doi.org/10.1016/j.neuroimage.2016.02.039

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
% Make paths to nii_xps folders
paths.nii_xps = cell(0,0);
for ndata = 1:numel(datasets_paths)
    datasets_path = datasets_paths{ndata};
    mdd_path = fullfile(datasets_path,'mdd');
    nifti_names = mdm_bruker_dir2expnams(mdd_path);
    
    for nnifti = 1:numel(nifti_names)
        nifti_name = nifti_names{nnifti};
        nifti_path = fullfile(mdd_path, nifti_name); 
        paths.nii_xps{1+numel(paths.nii_xps)} = fullfile(nifti_path,'nii_xps');
    end
end

% Prepare options
opt = mdm_opt();
opt = mio_opt(opt);
opt.verbose       = 1;
opt.do_overwrite  = 1;
opt.do_mask = 1; 
opt.mask.do_overwrite = 1;
opt.mask.threshold = .02;
opt.do_data2fit = 1; 
opt.do_fit2param = 1;
opt.mio.do_parfor = 1;
opt.dtd_covariance.do_clamping = 1;

methods = {'dtd_covariance'};
opt.dtd_covariance.fig_maps = {'s0','MD','meanshape','varsize'};
nii_name = 'data';

% Loop over datasets
for ndata = 1:numel(paths.nii_xps)

    for nmethod = numel(methods)
        method = methods{nmethod};

        % Connect to data
        clear s
        s.nii_fn = fullfile(paths.nii_xps{ndata}, [nii_name opt.nii_ext]);
        s.mask_fn = fullfile(paths.nii_xps{ndata}, [nii_name '_mask' opt.nii_ext]);
        s.xps = mdm_xps_load(mdm_xps_fn_from_nii_fn(s.nii_fn));

        % OUTPUT: define paths for data, fit parameters, and maps
        paths.mfs_fn   = fullfile(fileparts(paths.nii_xps{ndata}), method, 'mfs.mat');
        paths.dps_fn   = fullfile(fileparts(paths.nii_xps{ndata}), method, 'dps.mat');
        paths.nii_path   = fullfile(fileparts(paths.nii_xps{ndata}), method, 'nii');

        opt.mask.b0_ind = find(abs(s.xps.b-min(s.xps.b)) < .11e9);    

        % Run analysis
        pipename = [method '_pipe'];
        nii_fn = feval(pipename,s, paths, opt);
    end   
    
end

