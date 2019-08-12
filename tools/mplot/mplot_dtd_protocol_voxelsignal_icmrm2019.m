function mplot_dtd_protocol_voxelsignal(voxels, paths, opt)
% function mplot_dtd_protocol_voxelsignal(dps_fn, pdf_path, opt)

pix_i = voxels.i;
pix_j = voxels.j;
pix_k = voxels.k;

% pixel color
pix_col = [1 .3 .3
         0 0.8 0
         .3 .3 1];

% Font sizes etc
figscale = 1;
figwidth = figscale*30;
fs = figscale*24;
lw = figscale*3;
aspect = 3/4;
ms_max = figscale*30;
ms = figscale*6;

% Read data
[I,h] = mdm_nii_read(paths.nii_fn);
xps = mdm_xps_load(paths.xps_fn);
mfs = mdm_mfs_load(paths.mfs_fn);
dps = mdm_dps_load(paths.dps_fn);

%Merge solutions
dpar = [];
dperp = [];
theta = [];
phi = [];
%r2 = [];
w = [];
nn = 0;

bs_path = paths.bs_path;
bsno = msf_getdirno(bs_path);
%bsno = bsno(1:2);
for nbs = bsno
    mfs_fn   = fullfile(bs_path, num2str(nbs), 'mfs.mat');
    if exist(mfs_fn,'file')==2
        nbs;

        mfs = mdm_mfs_load(mfs_fn);

        m = mfs.m;
        sz = size(m);
        n = m(:,:,:,1);
        nn_temp = (sz(4)-1)/5;

        ind = false(sz(4),1);
        ind(2:5:end) = 1;

        dpar_temp = m(:,:,:,circshift(ind,0,1));
        dperp_temp = m(:,:,:,circshift(ind,1,1));
        theta_temp = m(:,:,:,circshift(ind,2,1));
        phi_temp = m(:,:,:,circshift(ind,3,1));
        w_temp = m(:,:,:,circshift(ind,4,1));


        dpar = cat(4,dpar,dpar_temp);
        dperp = cat(4,dperp,dperp_temp);
        theta = cat(4,theta,theta_temp);
        phi = cat(4,phi,phi_temp);
        w = cat(4,w,w_temp);
        nn = nn + nn_temp;
    end
end

dpar_4d = dpar;
dperp_4d = dperp;
theta_4d = theta;
phi_4d = phi;
w_4d = w;

%% Figure 3
% S0 with labeled pixels
figure(3), clf

nk = pix_k(1);

im2d_S0 = dps.s0(:,:,nk);
s0max = max(reshape(dps.s0,numel(dps.s0),1));

pixaspect = mfs.nii_h.pixdim(3)/mfs.nii_h.pixdim(2);

height = 1;
width = height/aspect/pixaspect;
%axh_S0_pixel = axes('position',[1-width 1-height-.12 width height]);
axh_S0_pixel = axes('position',[0 0 1 1]);
imagesc(im2d_S0'), hold on
colormap('gray')

set(axh_S0_pixel,'YDir','normal')
axis(axh_S0_pixel,'tight','off')

for axh_n = 1:numel(axh_S0_pixel)
    for pix_n = 1:numel(pix_i)
        h = plot(axh_S0_pixel(axh_n),pix_i(pix_n),pix_j(pix_n),'o');
        set(h,'Color',pix_col(pix_n,:),'LineWidth',1*lw,'MarkerSize',ms)
    end
end
datacursormode on

set(gcf, 'PaperUnits','centimeters','PaperPosition', 10*[0 0 1/pixaspect 1],'PaperSize', 10*[1/pixaspect 1]);
fig_fn = fullfile(paths.figs,['S0_LabeledVoxels']);
msf_mkdir(paths.figs);
eval(['print ' fig_fn ' -dpdf -loose'])

% Acquisition protocol
figure(4), clf
height = .1;
width = .4;
left = .12;
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
w_threshold = .1;
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
axh_dist_v = [];
ph_sign_v = [];
for pix_n = 1:numel(pix_i)
    ni = pix_i(pix_n);
    nj = pix_j(pix_n);
    left = left_v(pix_n);
    bottom = bottom_v(pix_n);

    s0 = im2d_S0(ni,nj,nk);

    signal = squeeze(I(ni,nj,nk,:))/s0;

    dpar = squeeze(dpar_4d(ni,nj,nk,:));
    dperp = squeeze(dperp_4d(ni,nj,nk,:));
    theta = squeeze(theta_4d(ni,nj,nk,:));
    phi = squeeze(phi_4d(ni,nj,nk,:));
    w = squeeze(w_4d(ni,nj,nk,:))/numel(bsno);
    s0 = sum(w);
    dtd = dtd_par2dist(dpar,dperp,theta,phi,w);
    dtd = dtd_sort(dtd);
    m = dtd;
    
    signal_fit = dtd_1d_fit2data(m, xps)/s0;
    %signal_fit = dtd_1d_fit2data(squeeze(mfs.m(ni,nj,nk,:))', xps)/s0;
    rmsd = sqrt(sum((signal_fit-signal).^2)/numel(signal));
    
    axh_pix = axes('position',[left bottom width height]);
    h = plot(1:numel(signal),signal(ind_sort)','o');
    set(h,'Color',pix_col(pix_n,:),'LineWidth',.5*lw,'MarkerSize',ms)
    hold on
    ph_sign = plot(1:numel(signal),signal_fit(ind_sort)','k.','MarkerSize',1.5*ms);
    ph_sign_v = [ph_sign_v; ph_sign];
    axh_pix_v = [axh_pix_v; axh_pix];

    %m = squeeze(mfs.m(ni,nj,nk,:))';
    %s0 = dps.s0(ni,nj,nk);
%     dtd = dtd_m2dtd(m);
    [n,dpar,dperp,theta,phi,w] = dtd_dist2par(dtd);
    if s0>0
        xcos = cos(phi).*sin(theta);
        ycos = sin(phi).*sin(theta);
        zcos = cos(theta);

        iso = tm_eigvals2iso([dpar dperp dperp]);
        fa = tm_eigvals2fa([dpar dperp dperp]);

        c.x = log10(iso);
        c.y = log10(dpar./dperp);
        c.ms = ms_max*sqrt(w*numel(bsno)/s0max);
        c.bright = fa;
        c.r = abs(xcos);
        c.g = abs(ycos);
        c.b = abs(zcos);
        
        [c.x c.y c.ms c.bright c.r c.g c.b];
        

        axh_dist = axes('position',[1-height/aspect-.06 bottom height/aspect height]);
        for nc = 1:n
            if w(nc) > w_threshold*s0/numel(bsno)
                h1 = plot(c.x(nc),c.y(nc),'o');
                hold on
                col = c.bright(nc)*[c.r(nc) c.g(nc) c.b(nc)];
                %set(h1,'MarkerSize',c.ms(nc),'Color',col,'MarkerFaceColor',col)
                set(h1,'MarkerSize',c.ms(nc),'Color',col,'MarkerFaceColor','none','LineWidth',.5*lw)
            end
        end
        axh_dist_v = [axh_dist_v; axh_dist];
    end
end

set(axh_b,'YLim',max(xps.b/1e9)*[-.1 1.1])
set(axh_bdelta,'YLim',[-.65 1.3],'YTick',[-.5 0 .5 1])
set(axh_theta,'YLim',180*[-.1 1.1],'YTick',[0:90:180])
set(axh_phi,'YLim',180*[-1.2 1.2],'YTick',[-180:90:180])
set([ph_prot_v; ph_sign_v],'Marker','.','LineStyle','none','Color','k','MarkerSize',1.5*ms)

set([axh_prot_v; axh_pix_v],'Box','off','TickDir','out','XLim',xps.n*[-.01 1.01],...
    'XTick',[0:100:500],'XTickLabel',[],'TickLength',.01*[1 1],'LineWidth',lw,'FontSize',fs)
set(axh_pix,'XTickLabel',[0:100:500])
set(axh_pix_v,'YLim',[-.05 1.05],'YTick',0:.5:1)

xmin = -9.5; xmax = -8; ymin =-1.5; ymax = 1.5;
set(axh_dist_v,'XLim',[xmin xmax]+.01*(xmax-xmin)*[-1 1], 'YLim',[ymin ymax]+.01*(ymax-ymin)*[-1 1],'XTick',-9.5:.5:-8,'YTick',-2:2,...
    'TickDir','out','TickLength',.025*[1 1],'LineWidth',lw,'FontSize',fs,'YAxisLocation','left','XTickLabel',[])
set(axh_dist,'XTickLabel',-9.5:.5:-8)

set(gcf, 'PaperUnits','centimeters','PaperPosition', figwidth*[0 0 1 1/aspect],'PaperSize', figwidth*[1 1/aspect]);
fig_fn = fullfile(paths.figs,['Figure3_protocol_pixelssignals']);
msf_mkdir(paths.figs);
eval(['print ' fig_fn ' -dpdf -loose'])



