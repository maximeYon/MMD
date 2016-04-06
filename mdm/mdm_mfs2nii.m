function mdm_mfs2nii(mfs_fn, o_path, figopt, opt)
% function mdm_mfs2nii(mfs_fn, o_path, figopt, opt)

if (nargin < 4), opt = []; end
 
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

for n = 1:numel(figopt.fig_maps);
    param = figopt.fig_maps{n};
    fn = fullfile(o_path,[figopt.fig_prefix '_' param opt.nii_ext]);
    mdm_nii_write(mfs.(param), fn, h);
end

