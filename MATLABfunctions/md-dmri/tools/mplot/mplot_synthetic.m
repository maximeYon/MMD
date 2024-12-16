function axh_v_tot = mplot_synthetic(method, dps, fig_fn, clim, opt)
% function mplot_synthetic(method, dps, fig_fn, clim, opt)

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
mask_nk = s0_nk>0;

% smax = max(dps.s0(:,:,nk),[],'all');
% smax = quantile(s0_nk(isfinite(s0_nk)),.99,'all');
smax = quantile(s0_nk(s0_nk>0),.999,'all');
clim_s0norm = clim.s0;
% clim.s0 = smax*clim.s0;
% clim.s2000 = max(dps.s2000(:,:,nk),[],'all')*clim.s2000;
% clim.s2000 = quantile(s2000_nk(isfinite(s0_nk)),.99,'all')*clim.s2000;

%Plot parameter maps
figure(2), clf

Nbins = min([numel(dps.bin) 3]);

if strcmp(method,'dtr2d')
%     plotfields.gray = {'s0';'s2000';'mdiso';'msddelta';'mr2';'vdiso';'vsddelta';'vr2'};
%     plotfields.hotcold = {'cvdisosddelta';'cvdisor2';'cvsddeltar2'};
%     plotfields.bin = {'mdiso';'msddelta';'mr2'};
%     plotfields.gray = {'s0';'s500';'s1000';'s2000';'s3000';'s4000';'s5000';'s8000'};
    plotfields.gray = {'s0';'s_b1000';'s_b2000';'s_b5000';'s_tr0020';'s_tr0050';'s_tr0100';'s_te0020';'s_te0050';'s_te0100'};
    plotfields.hotcold = {};
    plotfields.bin = {};
elseif strcmp(method,'dtr1d')
    plotfields.gray = {'s0';'s2000';'mdiso';'msddelta';'mr1';'vdiso';'vsddelta';'vr1'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvsddeltar1'};
    plotfields.bin = {'mdiso';'msddelta';'mr1'};
elseif strcmp(method,'dtr1r2d')
    plotfields.gray = {'s0';'s_b0500';'s_b1000';'s_b2000';'s_b0500_bdelta1';'s_b1000_bdelta1';'s_b2000_bdelta1';'s_tr1000';'s_tr0200';'s_tr0100';'s_te0020';'s_te0050';'s_te0100'};
    plotfields.hotcold = {};
    plotfields.bin = {};
elseif strcmp(method,'dtd')
    plotfields.gray = {'s0';'s500';'s1000';'s2000';'s3000';'s4000';'s5000';'s8000'};
%     plotfields.gray = {'s0';'s2000';'mdiso';'msddelta';'vdiso';'vsddelta'};
%     plotfields.gray = {'s0';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {};
    plotfields.bin = {};
elseif strcmp(method,'dtod')
    plotfields.gray = {'s0';'s_b1000';'s_b2000';'s_b5000'};
    plotfields.hotcold = {};
    plotfields.bin = {};
elseif strcmp(method,'dtor1r2d')
    plotfields.gray = {'s0';'s_b1000';'s_b2000';'s_b5000';'s_tr1000';'s_tr0100';'s_tr0050';'s_te0010';'s_te0020';'s_te0050'};
    plotfields.hotcold = {};
    plotfields.bin = {};
end


for nfield = 1:numel(plotfields.gray)
    field = plotfields.gray{nfield};
    dps_nk.(field) = dps.(field)(:,:,nk);
    clim.(field) = quantile(dps_nk.(field)(s0_nk>0),.999,'all')*clim_s0norm;    
end

pixaspect = dps.nii_h.pixdim(3)/dps.nii_h.pixdim(2);
imaspect = sz(2)/sz(1);

Nslices = numel(nk);
Nparams = numel(plotfields.gray) + numel(plotfields.hotcold) + Nbins*numel(plotfields.bin) + 0;
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
        im3d.bright = mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk);
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(im3d,position,[0 1]);
        axh_v_tot = [axh_v_tot; axh_v];
    end
end

% for nbin = 1:Nbins
%     clear im3d
%     im3d.r = dps.bin{nbin}.mdxx(:,:,nk);
%     im3d.g = dps.bin{nbin}.mdyy(:,:,nk);
%     im3d.b = dps.bin{nbin}.mdzz(:,:,nk);
%     im3d.bright = mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk);
%     position.left = position.left + dleft;
%     axh_v = mplot_slicescolumn(im3d,position,[0 1]);
%     axh_v_tot = [axh_v_tot; axh_v];
% end
% 
% clear im3d
% im3d.r = dps.bin{1}.f(:,:,nk);
% im3d.g = dps.bin{2}.f(:,:,nk);
% im3d.b = dps.bin{3}.f(:,:,nk);
% im3d.bright = mask(:,:,nk).*(dps.bin{1}.f(:,:,nk) + dps.bin{2}.f(:,:,nk) + dps.bin{3}.f(:,:,nk));
% clim = [0 1];
% position.left = position.left + dleft;
% axh_v = mplot_slicescolumn(im3d,position,clim);
% axh_v_tot = [axh_v_tot; axh_v];

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



