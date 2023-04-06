function axh_v_tot = mplot_technicolor(method, dps, fig_fn, clim, opt)
% function mplot_technicolor(method, dps, fig_fn, clim, opt)

if nargin < 5
    opt = mdm_opt();
end

sz = ones(1,3);
sz_temp = size(dps.s0);
sz(1:numel(sz_temp)) = sz_temp;

if isfield(opt,'k_range')
    nk = opt.k_range;
else
    nk = 1:sz(3); 
end

s0_nk = dps.s0(:,:,nk);
s_b0500_nk = dps.s_b0500(:,:,nk);
s_b1000_nk = dps.s_b1000(:,:,nk);
s_b2000_nk = dps.s_b2000(:,:,nk);
s_b3000_nk = dps.s_b3000(:,:,nk);
s_b5000_nk = dps.s_b5000(:,:,nk);
% smax = max(dps.s0(:,:,nk),[],'all');
% smax = quantile(s0_nk(isfinite(s0_nk)),.99,'all');

clim.s_b0500 = clim.s0;
clim.s_b1000 = clim.s0;
clim.s_b2000 = clim.s0;
clim.s_b3000 = clim.s0;
clim.s_b5000 = clim.s0;

smax = quantile(s0_nk(s0_nk>0),.999,'all');
clim.s0 = smax*clim.s0;
% clim.s2000 = max(dps.s2000(:,:,nk),[],'all')*clim.s2000;
% clim.s2000 = quantile(s2000_nk(isfinite(s0_nk)),.99,'all')*clim.s2000;
clim.s_b0500 = quantile(s_b0500_nk(s0_nk>0),.999,'all')*clim.s_b0500;
clim.s_b1000 = quantile(s_b1000_nk(s0_nk>0),.999,'all')*clim.s_b1000;
clim.s_b2000 = quantile(s_b2000_nk(s0_nk>0),.999,'all')*clim.s_b2000;
clim.s_b3000 = quantile(s_b3000_nk(s0_nk>0),.999,'all')*clim.s_b3000;
clim.s_b5000 = quantile(s_b5000_nk(s0_nk>0),.999,'all')*clim.s_b5000;

%Plot parameter maps
figure(2), clf

Nbins = min([numel(dps.bin) 3]);

if strcmp(method,'dtr2d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr2';'vdiso';'vsddelta';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor2';'cvsddeltar2'};
    plotfields.bin = {'mdiso';'msddelta';'mr2'};
elseif strcmp(method,'dtr1d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr1';'vdiso';'vsddelta';'vr1'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvsddeltar1'};
    plotfields.bin = {'mdiso';'msddelta';'mr1'};
elseif strcmp(method,'dtr1r2d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr1';'mr2';'vdiso';'vsddelta';'vr1';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2'};
    plotfields.bin = {'mdiso';'msddelta';'mr1';'mr2'};
elseif strcmp(method,'dtd')
    plotfields.gray = {'s0';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {'cvdisosddelta'};
    plotfields.bin = {'mdiso';'msddelta'};
elseif strcmp(method,'dtod')
    plotfields.gray = {'s0';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {'cvdisosddelta';'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
    plotfields.bin = {'mdiso';'msddelta'};
elseif strcmp(method,'dtor1r2d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr1';'mr2';'vdiso';'vsddelta';'vr1';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2';'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
    plotfields.bin = {'mdiso';'msddelta';'mr1';'mr2'};
end

pixaspect = dps.nii_h.pixdim(3)/dps.nii_h.pixdim(2);
imaspect = sz(2)/sz(1);

Nslices = numel(nk);
Nparams = numel(plotfields.gray) + numel(plotfields.hotcold) + Nbins*numel(plotfields.bin) + 4;
papersize = 3*[Nparams Nslices*pixaspect*imaspect];
position.dbottom = 1/Nslices;
dleft = 1/Nparams;

position.height = 1.01*position.dbottom; 
position.width = 1.01*dleft; 
position.left = -dleft;

mask = dps.s0 > clim.mask_threshold*smax;

axh_v_tot = [];

for c = 1:numel(plotfields.gray)
    im3d = double(dps.(plotfields.gray{c}));
    im3d(im3d>.99*max(clim.(plotfields.gray{c}))) = .99*max(clim.(plotfields.gray{c})); %Factor .99 to avoid black speckles in pdf        
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.gray{c}));
    for n = 1:numel(axh_v), colormap(axh_v(n),gray(256)), end
    axh_v_tot = [axh_v_tot; axh_v];
end

for c = 1:numel(plotfields.hotcold)
    im3d = double(dps.( plotfields.hotcold{c}));
    im3d(im3d>.99*max(clim.(plotfields.hotcold{c}))) = .99*max(clim.(plotfields.hotcold{c})); %Factor .99 to avoid black speckles in pdf
    im3d(im3d<.99*min(clim.(plotfields.hotcold{c}))) = .99*min(clim.(plotfields.hotcold{c})); %Factor .99 to avoid black speckles in pdf
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.hotcold{c}));
    for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(128)), end
    axh_v_tot = [axh_v_tot; axh_v];
end

for c = 1:numel(plotfields.bin)
    for nbin = 1:Nbins
        clear im3d
        cind = (dps.bin{nbin}.(plotfields.bin{c})(:,:,nk)-min(clim.(plotfields.bin{c})))...
            /(max(clim.(plotfields.bin{c}))-min(clim.(plotfields.bin{c})));
        im3d = dist_cind2rgb_jet(cind);
        im3d.bright = .99*mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk); %factor .99 added to avoid black voxels in pdf
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(im3d,position,[0 1]);
        axh_v_tot = [axh_v_tot; axh_v];
    end
end

for nbin = 1:Nbins
    clear im3d
    im3d.r = dps.bin{nbin}.mdxx(:,:,nk);
    im3d.g = dps.bin{nbin}.mdyy(:,:,nk);
    im3d.b = dps.bin{nbin}.mdzz(:,:,nk);
    im3d.bright = mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk);
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,[0 1]);
    axh_v_tot = [axh_v_tot; axh_v];
end

clear im3d
im3d.r = dps.bin{1}.f(:,:,nk);
im3d.g = dps.bin{2}.f(:,:,nk);
im3d.b = dps.bin{3}.f(:,:,nk);
im3d.bright = mask(:,:,nk).*(dps.bin{1}.f(:,:,nk) + dps.bin{2}.f(:,:,nk) + dps.bin{3}.f(:,:,nk));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
axh_v_tot = [axh_v_tot; axh_v];

if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
    set(axh_v_tot,'YDir','reverse')
elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
    set(axh_v_tot,'XDir','reverse','YDir','normal')
else
    set(axh_v_tot,'YDir','normal')
end

if isfield(opt,'roiplot')
    if isfield(opt.roiplot,'case_name')
        th_case = text(axh_v_tot(1),1,sz(2),opt.roiplot.case_name);
        set(th_case,'Color',[1 1 .999],'VerticalAlignment','top','FontSize',8)
    end
end

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

if ~isempty(fig_fn)
    msf_mkdir(fileparts(fig_fn));
    print(fig_fn,'-loose','-dpdf')
end



