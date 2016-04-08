function fn = mdm_param2nii(dps_fn, o_path, fig_opt, opt)
% function fn = mdm_param2nii(dps_fn, o_path, fig_opt, opt)

if (nargin < 4), opt = []; end

opt = mdm_opt(opt);

if (~opt.do_param2nii), disp('Returning (opt.do_param2nii == 0)'); end

fig_opt = msf_ensure_field(fig_opt, 'fig_cmaps', {}); % cmaps not required

dps = mdm_dps_load(dps_fn);
h   = dps.nii_h;
fn  = {};

% grayscale maps
for n = 1:numel(fig_opt.fig_maps);
    param = fig_opt.fig_maps{n};
    fn{end+1} = fullfile(o_path,[fig_opt.fig_prefix '_' param opt.nii_ext]);
    
    try
        mdm_nii_write(dps.(param), fn{end}, h);
    catch me
        disp(me.message);
    end
end

% color maps
is_colour = 1;

for n = 1:numel(fig_opt.fig_cmaps);
    param = fig_opt.fig_cmaps{n};
    col = fig_opt.fig_ccol{n};
    colnorm = fig_opt.fig_ccolnorm{n};
    
    % 2do: change to dps.(param)
    eval(['c.bright = dps.' param ';'])
    eval(['c.r = squeeze(abs(dps.' col '(:,:,:,1)))./dps.' colnorm ';'])
    eval(['c.g = squeeze(abs(dps.' col '(:,:,:,2)))./dps.' colnorm ';'])
    eval(['c.b = squeeze(abs(dps.' col '(:,:,:,3)))./dps.' colnorm ';'])
    
    Icol = zeros(3,size(c.bright,1),size(c.bright,2),size(c.bright,3));
    Icol(1,:,:,:) = c.bright.*c.r;
    Icol(2,:,:,:) = c.bright.*c.g;
    Icol(3,:,:,:) = c.bright.*c.b;
    Icol(isnan(Icol)) = 0;
    Icol(isinf(Icol)) = 0;
    Icol(Icol>1) = 1;
    Icol(Icol<0) = 0;
    
    fn{end+1} = fullfile(o_path,[fig_opt.fig_prefix '_' param '_' col '_rgb' opt.nii_ext]);
    mdm_nii_write(256*Icol, fn{end}, h, is_colour);
end

