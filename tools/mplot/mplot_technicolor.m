function mplot_technicolor(method, dps, fig_fn, clim)
% function mplot_technicolor(method, dps, fig_fn, clim)

%Plot parameter maps
figure(1), clf

smax = max(dps.s0(:));
clim.s0 = smax*clim.s0;
clim.s2000 = max(dps.s2000(:))*clim.s2000;

if strcmp(method,'dtr2d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr2';'vdiso';'vsddelta';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor2';'cvsddeltar2'};
    plotfields.bin = {'mdiso';'msddelta';'mr2'};
elseif strcmp(method,'dtr1d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr1';'vdiso';'vsddelta';'vr1'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvsddeltar1'};
    plotfields.bin = {'mdiso';'msddelta';'mr1'};
elseif strcmp(method,'dtd')
    plotfields.gray = {'s0';'s2000';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {'cvdisosddelta'};
    plotfields.bin = {'mdiso';'msddelta'};
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

axh_v_tot = [];

for c = 1:numel(plotfields.gray)
    im3d = double(dps.( plotfields.gray{c}));
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.gray{c}));
    for n = 1:numel(axh_v), colormap(axh_v(n),gray(256)), end
    axh_v_tot = [axh_v_tot; axh_v];
end

for c = 1:numel(plotfields.hotcold)
    im3d = double(dps.( plotfields.hotcold{c}));
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.hotcold{c}));
    for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(128)), end
    axh_v_tot = [axh_v_tot; axh_v];
end

for c = 1:numel(plotfields.bin)
    for nbin = 1:numel(dps.bin)
        clear im3d
        cind = (dps.bin{nbin}.(plotfields.bin{c})(:,:,nk)-min(clim.(plotfields.bin{c})))...
            /(max(clim.(plotfields.bin{c}))-min(clim.(plotfields.bin{c})));
        im3d = dist_cind2rgb_jet(cind);
        im3d.bright = mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk);
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(im3d,position,[0 1]);
        axh_v_tot = [axh_v_tot; axh_v];
    end
end

for nbin = 1:numel(dps.bin)
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
else
    set(axh_v_tot,'YDir','normal')
end

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

if ~isempty(fig_fn)
    msf_mkdir(fileparts(fig_fn));
    print(fig_fn,'-loose','-dpdf')
end



