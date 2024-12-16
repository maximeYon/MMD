function mplot_technicolor_nii(method, dps, o_path, clim, opt)
% function mplot_technicolor_nii(method, dps, o_path, clim, opt)

%Plot parameter maps

if nargin < 5
    opt = mdm_opt();
end

msf_mkdir(o_path)

sz = ones(1,3);
sz_temp = size(dps.s0);
sz(1:numel(sz_temp)) = sz_temp;

if isfield(opt,'k_range')
    nk = opt.k_range;
else
    nk = 1:sz(3); 
end

s0_nk = dps.s0(:,:,nk);
smax = quantile(s0_nk(s0_nk>0),.999,'all');
clim.s0 = smax*clim.s0;

%dtd
% plotfields.gray = {'s0';'s800';'s2000';'adc800';'mdiso';'msddelta';'vdiso';'vsddelta';'nmsdaniso';'nvdiso';'MD';'FA'};
plotfields.gray = {'s0';'s_b1000';'s_b2000';'s_b5000';'mdiso';'msddelta';'vdiso';'vsddelta';'nmsdaniso';'nvdiso';'MD';'FA'};
plotfields.hotcold = {'cvdisosddelta'};
plotfields.bin = {'mdiso';'msddelta'};
if strcmp(method,'dtr2d')
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b0200_te000';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b0500_te000';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te000';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te001';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te010';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te020';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te025';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te030';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te035';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te040';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te045';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te050';    
    plotfields.gray{1+numel(plotfields.gray)} = 'adc_b1000_te100';    
    plotfields.gray{1+numel(plotfields.gray)} = 'mr2';
    plotfields.gray{1+numel(plotfields.gray)} = 'vr2';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvdisor2';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvsddeltar2';
    plotfields.bin{1+numel(plotfields.bin)} = 'mr2';
elseif strcmp(method,'dtr1d')
    plotfields.gray{1+numel(plotfields.gray)} = 'mr1';
    plotfields.gray{1+numel(plotfields.gray)} = 'vr1';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvdisor1';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvsddeltar1';
    plotfields.bin{1+numel(plotfields.bin)} = 'mr1';
elseif strcmp(method,'dtr1r2d')
    plotfields.gray{1+numel(plotfields.gray)} = 's_te0010';
    plotfields.gray{1+numel(plotfields.gray)} = 's_te0020';
    plotfields.gray{1+numel(plotfields.gray)} = 's_te0050';
    plotfields.gray{1+numel(plotfields.gray)} = 's_te0100';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0020';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0050';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0100';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0200';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0500';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr1000';
    plotfields.gray{1+numel(plotfields.gray)} = 'mr1';
    plotfields.gray{1+numel(plotfields.gray)} = 'vr1';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvdisor1';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvsddeltar1';
    plotfields.bin{1+numel(plotfields.bin)} = 'mr1';
    plotfields.gray{1+numel(plotfields.gray)} = 'mr2';
    plotfields.gray{1+numel(plotfields.gray)} = 'vr2';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvdisor2';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvsddeltar2';
    plotfields.bin{1+numel(plotfields.bin)} = 'mr2';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvr1r2';
elseif strcmp(method,'dtod')
    plotfields.gray{1+numel(plotfields.gray)} = 'dmdisodnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'dmsddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'dvdisodnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'dvsddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'dcvdisosddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndmdisodnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndmsddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndvdisodnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndvsddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndcvdisosddeltadnu';
    plotfields.bin{1+numel(plotfields.bin)}   = 'dmdisodnu';
    plotfields.bin{1+numel(plotfields.bin)}   = 'dmsddeltadnu';
elseif strcmp(method,'dtor1r2d')
    plotfields.gray{1+numel(plotfields.gray)} = 's_te0010';
    plotfields.gray{1+numel(plotfields.gray)} = 's_te0020';
    plotfields.gray{1+numel(plotfields.gray)} = 's_te0050';
    plotfields.gray{1+numel(plotfields.gray)} = 's_te0100';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0020';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0050';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0100';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0200';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr0500';
    plotfields.gray{1+numel(plotfields.gray)} = 's_tr1000';
    plotfields.gray{1+numel(plotfields.gray)} = 'dmdisodnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'dmsddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'dvdisodnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'dvsddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'dcvdisosddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndmdisodnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndmsddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndvdisodnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndvsddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'ndcvdisosddeltadnu';
    plotfields.gray{1+numel(plotfields.gray)} = 'mr1';
    plotfields.gray{1+numel(plotfields.gray)} = 'vr1';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvdisor1';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvsddeltar1';
    plotfields.bin{1+numel(plotfields.bin)} = 'mr1';
    plotfields.gray{1+numel(plotfields.gray)} = 'mr2';
    plotfields.gray{1+numel(plotfields.gray)} = 'vr2';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvdisor2';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvsddeltar2';
    plotfields.bin{1+numel(plotfields.bin)} = 'mr2';
    plotfields.hotcold{1+numel(plotfields.hotcold)} = 'cvr1r2';
    plotfields.bin{1+numel(plotfields.bin)} = 'dmdisodnu';
    plotfields.bin{1+numel(plotfields.bin)} = 'dmsddeltadnu';
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
    if isfield(dps,plotfields.gray{c})
    tmp_fn = fullfile(o_path, [method '_' plotfields.gray{c} opt.nii_ext]); 
    im3d = mask.*dps.(plotfields.gray{c});
    nii_fn = mdm_nii_write(im3d, tmp_fn, dps.nii_h, 0);
    end
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
    
    
    tmp_fn = fullfile(o_path, [method '_' plotfields.hotcold{c} '_gray' opt.nii_ext]);
    im3d = mask.*dps.(plotfields.hotcold{c});
    nii_fn = mdm_nii_write(im3d, tmp_fn, dps.nii_h, 0);
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

        tmp_fn = fullfile(o_path, [method '_' plotfields.bin{c} '_bin' num2str(nbin) '_gray' opt.nii_ext]);
        
        clear im3d
        im3d = mask.*dps.bin{nbin}.(plotfields.bin{c});
        nii_fn = mdm_nii_write(im3d, tmp_fn, dps.nii_h, 0);        
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

% im3d.bright = mask.*(dps.bin{1}.f + dps.bin{2}.f + dps.bin{3}.f);
im3d.bright = mask;

I = zeros([3 msf_size(im3d.r,3)]);
I(1,:,:,:) = im3d.bright .* im3d.r ./c_norm;
I(2,:,:,:) = im3d.bright .* im3d.g ./c_norm;
I(3,:,:,:) = im3d.bright .* im3d.b ./c_norm;

I(isnan(I)) = 0;
I(isinf(I)) = 0;

I( I(:) > 1 ) = 1;
I( I(:) < 0 ) = 0;

mdm_nii_write(255*I, tmp_fn, dps.nii_h, 1);


tmp_fn = fullfile(o_path, [method '_s0_fractions' opt.nii_ext]);

clear im3d
im3d.r = dps.bin{1}.f;
im3d.g = dps.bin{2}.f;
im3d.b = dps.bin{3}.f;
c_norm = zeros([3 msf_size(im3d.r,3)]);
c_norm(1,:,:,:) = im3d.r;
c_norm(2,:,:,:) = im3d.g;
c_norm(3,:,:,:) = im3d.b;
c_norm = squeeze(max(c_norm,[],1));

im3d.bright = mask.*dps.s0/smax;

I = zeros([3 msf_size(im3d.r,3)]);
I(1,:,:,:) = im3d.bright .* im3d.r ./c_norm;
I(2,:,:,:) = im3d.bright .* im3d.g ./c_norm;
I(3,:,:,:) = im3d.bright .* im3d.b ./c_norm;

I(isnan(I)) = 0;
I(isinf(I)) = 0;

I( I(:) > 1 ) = 1;
I( I(:) < 0 ) = 0;

mdm_nii_write(255*I, tmp_fn, dps.nii_h, 1);


% Standard FA DEC
tmp_fn = fullfile(o_path, [method '_FA_u_rgb' opt.nii_ext]);

clear im3d
im3d.r = abs(dps.u(:,:,:,1));
im3d.g = abs(dps.u(:,:,:,2));
im3d.b = abs(dps.u(:,:,:,3));
c_norm = mask;

im3d.bright = dps.FA;

I = zeros([3 msf_size(im3d.r,3)]);
I(1,:,:,:) = im3d.bright .* im3d.r ./c_norm;
I(2,:,:,:) = im3d.bright .* im3d.g ./c_norm;
I(3,:,:,:) = im3d.bright .* im3d.b ./c_norm;

I(isnan(I)) = 0;
I(isinf(I)) = 0;

I( I(:) > 1 ) = 1;
I( I(:) < 0 ) = 0;

mdm_nii_write(255*I, tmp_fn, dps.nii_h, 1);

%% Fraction and md xx yy zz as double bin1, bin2, bin3
for nbin = 1:numel(dps.bin)
tmp_fn = fullfile(o_path, [method '_f_bin' num2str(nbin) '_gray' opt.nii_ext]);
nii_fn = mdm_nii_write(dps.bin{1, nbin}.f  , tmp_fn, dps.nii_h, 0);  

tmp_fn = fullfile(o_path, [method '_mdxx_bin' num2str(nbin) '_gray' opt.nii_ext]);
nii_fn = mdm_nii_write(dps.bin{nbin}.mdxx  , tmp_fn, dps.nii_h, 0);  
tmp_fn = fullfile(o_path, [method '_mdyy_bin' num2str(nbin) '_gray' opt.nii_ext]);
nii_fn = mdm_nii_write(dps.bin{nbin}.mdyy  , tmp_fn, dps.nii_h, 0); 
tmp_fn = fullfile(o_path, [method '_mdzz_bin' num2str(nbin) '_gray' opt.nii_ext]);
nii_fn = mdm_nii_write(dps.bin{nbin}.mdzz  , tmp_fn, dps.nii_h, 0); 
end

