function mplot_technicolor_nii(method, dps, o_path, clim, opt)
% function mplot_technicolor_nii(method, dps, o_path, clim, opt)

%Plot parameter maps

smax = max(dps.s0(:));
clim.s0 = smax*clim.s0;

%dtd
plotfields.gray = {'s0';'s1000';'s2000';'s3000';'mdiso';'msddelta';'vdiso';'vsddelta';'nmsdaniso';'nvdiso'};
plotfields.hotcold = {'cvdisosddelta'};
plotfields.bin = {'mdiso';'msddelta'};
if strcmp(method,'dtr2d')
    plotfields.gray = {plotfields.gray;'mr2';'vr2'};
    plotfields.hotcold = {plotfields.hotcold;'cvdisor2';'cvsddeltar2'};
    plotfields.bin = {plotfields.bin;'mr2'};
elseif strcmp(method,'dtr1d')
    plotfields.gray = {plotfields.gray;'mr1';'vr1'};
    plotfields.hotcold = {plotfields.hotcold;'cvdisor1';'cvsddeltar1'};
    plotfields.bin = {plotfields.bin;'mr1'};
end

sz = ones(1,3);
sz_temp = size(dps.s0);
sz(1:numel(sz_temp)) = sz_temp;
pixaspect = dps.nii_h.pixdim(3)/dps.nii_h.pixdim(2);
imaspect = sz(2)/sz(1);

nk = 1:sz(3); %nk = [2 6 11];
Nslices = numel(nk);
Nparams = numel(plotfields.gray) + numel(plotfields.hotcold) + numel(dps.bin)*numel(plotfields.bin) + 4;
papersize = 3*[Nparams Nslices*pixaspect*imaspect];
position.dbottom = 1/Nslices;
dleft = 1/Nparams;

position.height = 1.01*position.dbottom; 
position.width = 1.01*dleft; 
position.left = -dleft;

mask = dps.s0 > clim.mask_threshold*smax;

for c = 1:numel(plotfields.gray)
    tmp_fn = fullfile(o_path, [method '_' plotfields.gray{c} opt.nii_ext]);
    im3d = mask.*dps.(plotfields.gray{c});
    nii_fn = mdm_nii_write(im3d, tmp_fn, dps.nii_h, 0);
end

for c = 1:numel(plotfields.hotcold)
    tmp_fn = fullfile(o_path, [method '_' plotfields.hotcold{c} opt.nii_ext]);
    
    cmap = mplot_cmaphotcold(64);
    
    clim_temp = clim.(plotfields.hotcold{c});
    im3d = mask.*dps.(plotfields.hotcold{c});
    im3d(im3d<clim_temp(1)) = clim_temp(1);
    im3d(im3d>clim_temp(2)) = clim_temp(2);
    im3d_ind = ceil((im3d-clim_temp(1))/(clim_temp(2)-clim_temp(1))*size(cmap,1));
        
    I = zeros([3 msf_size(im3d,3)]);
    for k = nk
        RGB = ind2rgb(im3d_ind(:,:,k),cmap);
       
        I(1,:,:,k) = RGB(:,:,1);
        I(2,:,:,k) = RGB(:,:,2);
        I(3,:,:,k) = RGB(:,:,3);
    end

    I(isnan(I)) = 0;
    I(isinf(I)) = 0;

    I( I(:) > 1 ) = 1;
    I( I(:) < 0 ) = 0;

    mdm_nii_write(255*I, tmp_fn, dps.nii_h, 1);    
end

for c = 1:numel(plotfields.bin)
    for nbin = 1:numel(dps.bin)
        tmp_fn = fullfile(o_path, [method '_' plotfields.bin{c} '_bin' num2str(nbin) opt.nii_ext]);
        
        clear im3d
        cind = (dps.bin{nbin}.(plotfields.bin{c})(:,:,nk)-min(clim.(plotfields.bin{c})))...
            /(max(clim.(plotfields.bin{c}))-min(clim.(plotfields.bin{c})));
        im3d = dist_cind2rgb_jet(cind);
        im3d.bright = mask.*dps.bin{nbin}.f;
        
        I = zeros([3 msf_size(im3d.r,3)]);
        I(1,:,:,:) = im3d.bright .* im3d.r;
        I(2,:,:,:) = im3d.bright .* im3d.g;
        I(3,:,:,:) = im3d.bright .* im3d.b;

        I(isnan(I)) = 0;
        I(isinf(I)) = 0;

        I( I(:) > 1 ) = 1;
        I( I(:) < 0 ) = 0;
    
        mdm_nii_write(255*I, tmp_fn, dps.nii_h, 1);
    end
end

for nbin = 1:numel(dps.bin)
    tmp_fn = fullfile(o_path, [method '_mdii_bin' num2str(nbin) opt.nii_ext]);
    clear im3d
    im3d.r = dps.bin{nbin}.mdxx;
    im3d.g = dps.bin{nbin}.mdyy;
    im3d.b = dps.bin{nbin}.mdzz;
    c_norm = zeros([3 msf_size(im3d.r,3)]);
    c_norm(1,:,:,:) = im3d.r;
    c_norm(2,:,:,:) = im3d.g;
    c_norm(3,:,:,:) = im3d.b;
    c_norm = squeeze(max(c_norm,[],1));
    
    im3d.bright = mask.*dps.bin{nbin}.f;
    
    I = zeros([3 msf_size(im3d.r,3)]);
    I(1,:,:,:) = im3d.bright .* im3d.r ./c_norm;
    I(2,:,:,:) = im3d.bright .* im3d.g ./c_norm;
    I(3,:,:,:) = im3d.bright .* im3d.b ./c_norm;

    I(isnan(I)) = 0;
    I(isinf(I)) = 0;

    I( I(:) > 1 ) = 1;
    I( I(:) < 0 ) = 0;

    mdm_nii_write(255*I, tmp_fn, dps.nii_h, 1);
end

tmp_fn = fullfile(o_path, [method '_fractions' opt.nii_ext]);

clear im3d
im3d.r = dps.bin{1}.f;
im3d.g = dps.bin{2}.f;
im3d.b = dps.bin{3}.f;
c_norm = zeros([3 msf_size(im3d.r,3)]);
c_norm(1,:,:,:) = im3d.r;
c_norm(2,:,:,:) = im3d.g;
c_norm(3,:,:,:) = im3d.b;
c_norm = squeeze(max(c_norm,[],1));

im3d.bright = mask.*(dps.bin{1}.f + dps.bin{2}.f + dps.bin{3}.f);

I = zeros([3 msf_size(im3d.r,3)]);
I(1,:,:,:) = im3d.bright .* im3d.r ./c_norm;
I(2,:,:,:) = im3d.bright .* im3d.g ./c_norm;
I(3,:,:,:) = im3d.bright .* im3d.b ./c_norm;

I(isnan(I)) = 0;
I(isinf(I)) = 0;

I( I(:) > 1 ) = 1;
I( I(:) < 0 ) = 0;

mdm_nii_write(255*I, tmp_fn, dps.nii_h, 1);



