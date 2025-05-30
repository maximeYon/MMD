function mplot_dtd_protocol_voxelsignal(voxels, paths, opt)
% function mplot_dtd_protocol_voxelsignal(dps_fn, pdf_path, opt)

pix_i = voxels.i;
pix_j = voxels.j;
pix_k = voxels.k;

% pixel color
pix_col = [0 0.8 0
         1 .3 .3
         0 0 1];

% Font sizes etc
figscale = 2;
figwidth = figscale*17.78;
fs = figscale*6;
lw = figscale*1;
aspect = 1.618;
ms_max = figscale*20;
ms = figscale*3;

% Read data
[I,h] = mdm_nii_read(paths.nii_fn);
xps = mdm_xps_load(paths.xps_fn);
mfs = mdm_mfs_load(paths.mfs_fn);
dps = mdm_dps_load(paths.dps_fn);

%% Figure 3
% S0 with labeled pixels
figure(3), clf

nk = pix_k(1);

im2d_S0 = dps.s0(:,:,nk);
s0max = max(reshape(dps.s0,numel(dps.s0),1));

height = .27;
width = height/aspect;
axh_S0_pixel = axes('position',[1-width 1-height-.12 width height]);
imagesc(im2d_S0'), hold on
colormap('gray')

set(axh_S0_pixel,'YDir','normal')
axis(axh_S0_pixel,'equal','tight','off')

for axh_n = 1:numel(axh_S0_pixel)
    for pix_n = 1:numel(pix_i)
        h = plot(axh_S0_pixel(axh_n),pix_i(pix_n),pix_j(pix_n),'o');
        set(h,'Color',pix_col(pix_n,:),'LineWidth',1*lw,'MarkerSize',ms)
    end
end


% Acquisition protocol
height = .1;
width = .7;
left = .05;
bottom_v = .52 + 1.14*[0*height 1*height 2*height 3*height];

xps_array = [round(xps.b,-7), round(xps.b_delta,1), xps.theta/pi*180, xps.phi/pi*180];
[xps_array_sort, ind_sort] = sortrows(xps_array,[1 2 3 4]);

axh_b = axes('position',[left bottom_v(4) width height]);
ph_b = plot(axh_b,1:xps.n,xps.b(ind_sort)/1e9);
axh_bdelta = axes('position',[left bottom_v(3) width height]);
ph_bdelta = plot(axh_bdelta,1:xps.n,xps.b_delta(ind_sort));
axh_theta = axes('position',[left bottom_v(2) width height]);
ph_theta = plot(axh_theta,1:xps.n,xps.theta(ind_sort)/pi*180);
axh_phi = axes('position',[left bottom_v(1) width height]);
ph_phi = plot(axh_phi,1:xps.n,xps.phi(ind_sort)/pi*180);

axh_prot_v = [axh_b; axh_bdelta; axh_theta; axh_phi];
ph_prot_v = [ph_b; ph_bdelta; ph_theta; ph_phi];

% Signals and DTDs for selected pixels
w_threshold = .01;
s0_threshold = .1;

height = .12;
left_v = left + [0 0 0];
bottom_v = .06 + 1.15*[2*height height 0];

dmin = opt.dtd.dmin;
dmax = opt.dtd.dmax;
ratiomin = dmin/dmax;
ratiomax = dmax/dmin;

xmin = log10(dmin);
xmax = log10(dmax);
ymin = log10(1/ratiomax);
ymax = log10(ratiomax);

axh_pix_v = [];
axh_logs_v = [];
axh_dist_v = [];
for pix_n = 1:numel(pix_i)
    ni = pix_i(pix_n);
    nj = pix_j(pix_n);
    left = left_v(pix_n);
    bottom = bottom_v(pix_n);

    s0 = im2d_S0(ni,nj,nk);

    signal = squeeze(I(ni,nj,nk,:))/s0;
    signal_fit = dtd_1d_fit2data(squeeze(mfs.m(ni,nj,nk,:))', xps)/s0;
    rmsd = sqrt(sum((signal_fit-signal).^2)/numel(signal));
    
    axh_pix = axes('position',[left bottom width height]);
    h = plot(1:numel(signal),signal(ind_sort)','o');
    set(h,'Color',pix_col(pix_n,:),'LineWidth',1*lw,'MarkerSize',ms)
    hold on
    plot(1:numel(signal),signal_fit(ind_sort)','k.','MarkerSize',1.2*ms)
    axh_pix_v = [axh_pix_v; axh_pix];
    
    axh_logs = axes('position',[1-height/aspect-.14 bottom height/aspect height]);
    h = plot(xps.b/1e9,log10(signal'),'o');
    set(h,'Color',pix_col(pix_n,:),'LineWidth',1*lw,'MarkerSize',ms)
    hold on
    plot(xps.b/1e9,log10(signal_fit'),'k.','MarkerSize',1.2*ms)
    axh_logs_v = [axh_logs_v; axh_logs];

    m = squeeze(mfs.m(ni,nj,nk,:))';
    s0 = dps.s0(ni,nj,nk);
    dtd = dtd_m2dtd(m);
    [n,par,perp,theta,phi,w] = dtd_dist2par(dtd);
    if n>0
        xcos = cos(phi).*sin(theta);
        ycos = sin(phi).*sin(theta);
        zcos = cos(theta);

        iso = tm_eigvals2iso([par perp perp]);
        fa = tm_eigvals2fa([par perp perp]);

        c.x = log10(iso);
        c.y = log10(par./perp);
        c.ms = ms_max*sqrt(w/s0max);
        c.bright = fa;
        c.r = abs(xcos);
        c.g = abs(ycos);
        c.b = abs(zcos);

        axh_dist = axes('position',[1-height/aspect-.04 bottom height/aspect height]);
        for nc = 1:n
            if w(nc) > w_threshold*s0
                h1 = plot(c.x(nc),c.y(nc),'o');
                hold on
                col = c.bright(nc)*[c.r(nc) c.g(nc) c.b(nc)];
                %set(h1,'MarkerSize',c.ms(nc),'Color',col,'MarkerFaceColor',col)
                set(h1,'MarkerSize',c.ms(nc),'Color',col,'MarkerFaceColor','none','LineWidth',.75*lw)
            end
        end
        axh_dist_v = [axh_dist_v; axh_dist];
    end
end

set(axh_b,'YLim',max(xps.b/1e9)*[-.1 1.1])
set(axh_bdelta,'YLim',[-.65 1.3],'YTick',[-.5 0 .5 1])
set(axh_theta,'YLim',180*[-.1 1.1],'YTick',[0:90:180])
set(axh_phi,'YLim',180*[-1.2 1.2],'YTick',[-180:90:180])
set(ph_prot_v,'Marker','.','LineStyle','none','Color','k','MarkerSize',1.2*ms)

set([axh_prot_v; axh_pix_v],'Box','off','TickDir','out','XLim',xps.n*[-.01 1.01],...
    'XTick',[0:100:500],'XTickLabel',[],'TickLength',.005*[1 1],'LineWidth',lw,'FontSize',fs)
set(axh_pix,'XTickLabel',[0:100:500])
set(axh_pix_v,'YLim',[-.05 1.05],'YTick',0:.5:1)

set([axh_logs_v],'Box','off','TickDir','out','XLim',max(xps.b/1e9)*[-.01 1.01],...
    'XTick',[0:2:10],'XTickLabel',[],'TickLength',.025*[1 1],'LineWidth',lw,'FontSize',fs)
set([axh_logs],'XTickLabel',[0:2:10])

set(axh_dist_v,'XLim',[xmin xmax]+.01*(xmax-xmin)*[-1 1], 'YLim',[ymin ymax]+.01*(ymax-ymin)*[-1 1],'XTick',-11:-8,'YTick',-2:2,...
    'TickDir','out','TickLength',.025*[1 1],'LineWidth',lw,'FontSize',fs,'YAxisLocation','right','XTickLabel',[])
set(axh_dist,'XTickLabel',-11:-8)

set(gcf, 'PaperUnits','centimeters','PaperPosition', figwidth*[0 0 1 1/aspect],'PaperSize', figwidth*[1 1/aspect]);
fig_fn = fullfile(paths.figs,['Figure3_protocol_pixelssignals']);
msf_mkdir(paths.figs);
eval(['print ' fig_fn ' -dpdf -loose'])



