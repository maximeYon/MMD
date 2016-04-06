function mdm_mfs2cnii(mfs_fn, o_path, figopt, opt)
% function mdm_mfs2cnii(mfs_fn, o_path, param_v, prefix, opt)

if (nargin < 4), opt = []; end
 
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

is_colour = 1;

for n = 1:numel(figopt.fig_cmaps);
    param = figopt.fig_cmaps{n};
    col = figopt.fig_ccol{n};
    colnorm = figopt.fig_ccolnorm{n};

    % 2do: change to mfs.(param)
    eval(['c.bright = mfs.' param ';'])
    eval(['c.r = squeeze(abs(mfs.' col '(:,:,:,1)))./mfs.' colnorm ';'])
    eval(['c.g = squeeze(abs(mfs.' col '(:,:,:,2)))./mfs.' colnorm ';'])
    eval(['c.b = squeeze(abs(mfs.' col '(:,:,:,3)))./mfs.' colnorm ';'])

    Icol = zeros(3,size(c.bright,1),size(c.bright,2),1);
    Icol(1,:,:,:) = c.bright.*c.r;
    Icol(2,:,:,:) = c.bright.*c.g;
    Icol(3,:,:,:) = c.bright.*c.b;
    Icol(isnan(Icol)) = 0;
    Icol(isinf(Icol)) = 0;
    Icol(Icol>1) = 1;
    Icol(Icol<0) = 0;

    fn = fullfile(o_path,[figopt.fig_prefix '_' param '_' col '_rgb' opt.nii_ext]);
    mdm_nii_write(256*Icol, fn, h, is_colour);
end

