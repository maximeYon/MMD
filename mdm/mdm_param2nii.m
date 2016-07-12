function fn = mdm_param2nii(dps_fn, o_path, fig_opt, opt)
% function fn = mdm_param2nii(dps_fn, o_path, fig_opt, opt)

if (nargin < 4), opt = []; end

opt = mdm_opt(opt);

if (opt.do_param2nii == 0)
    disp('Returning (opt.do_param2nii == 0)'); 
    return;
end

% init: no not require any figures to be made
fig_opt = msf_ensure_field(fig_opt, 'fig_maps', {}); 
fig_opt = msf_ensure_field(fig_opt, 'fig_cmaps', {}); 

dps = mdm_dps_load(dps_fn);
fn  = {};

% grayscale maps
for n = 1:numel(fig_opt.fig_maps);
    
    param = fig_opt.fig_maps{n};
    
    tmp_fn = fullfile(o_path, [fig_opt.fig_prefix '_' param opt.nii_ext]);
    
    try
        mdm_nii_write(dps.(param), tmp_fn, dps.nii_h);
    catch me
        disp(me.message);
    end
    
    fn{end+1} = tmp_fn;
end

% color maps
for n = 1:numel(fig_opt.fig_cmaps);
    
    param   = fig_opt.fig_cmaps{n};
    col     = fig_opt.fig_ccol{n};
    colnorm = fig_opt.fig_ccolnorm{n};
    
    c.bright = dps.(param);
    c.r = squeeze(abs(dps.(col)(:,:,:,1))) ./ dps.(colnorm);
    c.g = squeeze(abs(dps.(col)(:,:,:,2))) ./ dps.(colnorm);
    c.b = squeeze(abs(dps.(col)(:,:,:,3))) ./ dps.(colnorm);
    
    I = zeros([3 msf_size(c.r,3)]); 
    I(1,:,:,:) = c.bright .* c.r;
    I(2,:,:,:) = c.bright .* c.g;
    I(3,:,:,:) = c.bright .* c.b;
    
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    
    I( I(:) > 1 ) = 1;
    I( I(:) < 0 ) = 0;
    
    tmp_fn = fullfile(o_path, [fig_opt.fig_prefix '_' param '_' col '_rgb' opt.nii_ext]);
    
    try
        mdm_nii_write(255*I, tmp_fn, dps.nii_h, 1);
    catch me
        disp(me.message);
    end
    
    fn{end+1} = tmp_fn;
end

