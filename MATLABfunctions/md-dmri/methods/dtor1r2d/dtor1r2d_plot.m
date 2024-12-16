function m = dtor1r2d_plot(S, xps, axh, axh2, opt)
% function dtor1r2d_plot(S, xps, axh, axh2)

if (nargin < 4), axh2 = []; end
if (nargin < 5), opt = []; end

method = 'dtor1r2d';

opt = mdm_opt(opt);
opt = dtd_opt(opt);
opt = feval([method '_opt'],opt);
opt = mplot_opt(opt);

opt.(method).dmin = .02/max(xps.b);
% opt.(method).dmax = 1/min(xps.b(xps.b>0));
opt.(method).dmax = 5e-9;
opt.(method).rmin = .01*min(xps.momega(xps.momega>0)/2/pi);
opt.(method).rmax = 100*max(xps.momega/2/pi);
opt.(method).maps_omega = min(xps.momega(xps.b>.5*max(xps.b)));
opt.(method).r1min = .1/max(xps.tr(xps.tr<20));
opt.(method).r1max = 2/min(xps.tr);
opt.(method).r2min = .01/max(xps.te);
opt.(method).r2max = 2/min(xps.te);
    
ms = 5;
fs = 10;
lw = 1;

S = abs(S);

% Show signal and fit
m = feval('dtor1r2d_1d_data2fit', S, xps, opt);
%m = [2 1e-9*[1 1] 0 0 2 8 3e8 1e-10*[1 1] 0 0 2 8 3e8];
S_fit = feval('dtor1r2d_1d_fit2data', m, xps)';

% Get dps from dtor1r2d_4d_fit2param
mfs.m = zeros(1,1,1,numel(m)); mfs.m(1,1,1,:) = m;
dps = dtor1r2d_4d_fit2param(mfs.m);


dmin = opt.dtor1r2d.dmin;
dmax = opt.dtor1r2d.dmax;
ratiomin = dmin/dmax;
ratiomax = dmax/dmin;
r1min = opt.dtor1r2d.r1min
r1max = opt.dtor1r2d.r1max;
r2min = opt.dtor1r2d.r2min;
r2max = opt.dtor1r2d.r2max

xmin = log10(dmin);
xmax = log10(dmax);
ymin = log10(1/ratiomax);
ymax = log10(ratiomax);
zmin = log10(r1min);
zmax = log10(r1max);

opt.(method).maps_omega/2/pi
[n,dpar,dperp,theta,phi,r1,r2,w] = dtor1r2d_dist2par(dtor1r2d_m2dtor1r2d(m));
[dpar,dperp,theta,phi,r1,r2,w] = dtor1r2d_4d_m2parso(mfs.m,opt.(method).maps_omega);
% [dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(squeeze(dpar),squeeze(dperp),squeeze(theta),squeeze(phi));
% dtd_1x6 = [dxx dyy dzz sqrt(2)*[dxy dxz dyz]];

diso = (dpar+2*dperp)/3;
dratio = dpar./dperp;
s0 = sum(w);

cla(axh);
hold(axh, 'off');
%figure(2), clf, axh = axes('position',[.15 .15 .7 .7]);
h1 = plot(axh,1:xps.n,S,'ro',1:xps.n,S_fit,'k.');
set(h1,'MarkerSize',ms,'LineWidth',lw)
axis(axh,'tight')
set(axh,'XLim',xps.n*[-.1 1.1], 'YLim',s0*[-.1 1.1],...
    'Box','off','TickDir','out','TickLength',.02*[1 1],...
'FontSize',fs,'LineWidth',lw)
set(axh,'YLim',max([S_fit(:); S(:)])*[-.1 1.1])
xlabel(axh,'Acq number','FontSize',fs)
ylabel(axh,'Signal','FontSize',fs)
title(axh,['<\omega> / 2\pi = ' num2str(min(xps.momega(xps.b>.5*max(xps.b))/2/pi),3) ' - ' num2str(max(xps.momega(xps.b>.5*max(xps.b))/2/pi),3) ' Hz'])
%return
fh2 = figure(2);
clf(fh2)
fh2_axh = copyobj(axh,fh2);
set(fh2_axh,'position',[.05 .15 .9 .8])

cla(axh2);
hold(axh2, 'on');
axis(axh2,'on')

if n>0

    if (1)
        axh_r1 = axh;
        dist_d.x = log10(diso);
        dist_d.y = log10(dratio);
        dist_d.z = log10(r1);
        dist_d.a = 500*w/s0;
        dist_d.r = abs(cos(phi).*sin(theta));
        dist_d.g = abs(sin(phi).*sin(theta));
        dist_d.b = abs(cos(theta));
        dist_d.bright = (abs(squeeze(dpar)-squeeze(dperp))./max([squeeze(dpar) squeeze(dperp)],[],2)).^2;

        contourpars.Nx = 50; contourpars.Ny = contourpars.Nx; contourpars.Nz = contourpars.Nx; contourpars.Nlevels = 5;
        axpars.xmin = log10(dmin)-.2; axpars.xmax = log10(dmax)+.2; axpars.ymin = log10(ratiomin)-.2; axpars.ymax = log10(ratiomax)+.2; axpars.zmin = log10(r1min)-.2; axpars.zmax = log10(r1max)+.2; 

        [hax,hscatter,hcontour] = dist_3d_scattercontourplot(axh_r1,dist_d,contourpars,axpars);

        set(axh_r1,'LineWidth',lw,'FontSize',fs)
        set(axh_r1,'XTick',-12:.5:-8,'YTick',-3:1:3,'ZTick',-2:.5:2,'XGrid','on','YGrid','on','ZGrid','on','Projection','perspective')

        xlabel(axh_r1,{'log_{10}(\itD\rm_{iso} / m^2s^-^1)'; '"size"'},'FontSize',fs) 
        ylabel(axh_r1,{'log_{10}(\itD\rm_{A} / \itD\rm_{R})'; '"shape"'},'FontSize',fs);
        zlabel(axh_r1,'log_{10}(\it{R}\rm_1 / s^-^1)','FontSize',fs)
        title(axh_r1,['orientation, [RGB]=[xyz]'],'FontSize',fs,'FontWeight','normal');
    end
    
    if (1)
        dist_d.x = log10(diso);
        dist_d.y = log10(dratio);
        dist_d.z = log10(r2);
        dist_d.a = 500*w/s0;
        dist_d.r = abs(cos(phi).*sin(theta));
        dist_d.g = abs(sin(phi).*sin(theta));
        dist_d.b = abs(cos(theta));
        dist_d.bright = (abs(squeeze(dpar)-squeeze(dperp))./max([squeeze(dpar) squeeze(dperp)],[],2)).^2;

        contourpars.Nx = 50; contourpars.Ny = contourpars.Nx; contourpars.Nz = contourpars.Nx; contourpars.Nlevels = 5;
        axpars.xmin = log10(dmin)-.2; axpars.xmax = log10(dmax)+.2; axpars.ymin = log10(ratiomin)-.2; axpars.ymax = log10(ratiomax)+.2; axpars.zmin = log10(r2min)-.2; axpars.zmax = log10(r2max)+.5; 

        [hax,hscatter,hcontour] = dist_3d_scattercontourplot(axh2,dist_d,contourpars,axpars);

        set(axh2,'LineWidth',lw,'FontSize',fs)
        set(axh2,'XTick',-12:.5:-8,'YTick',-3:1:3,'ZTick',-1:.5:3,'XGrid','on','YGrid','on','ZGrid','on','Projection','perspective')

        xlabel(axh2,{'log_{10}(\itD\rm_{iso} / m^2s^-^1)'; '"size"'},'FontSize',fs) 
        ylabel(axh2,{'log_{10}(\itD\rm_{A} / \itD\rm_{R})'; '"shape"'},'FontSize',fs);
        zlabel(axh2,'log_{10}(\it{R}\rm_2 / s^-^1)','FontSize',fs)
        title(axh2,['orientation, [RGB]=[xyz]'],'FontSize',fs,'FontWeight','normal');
    end

end

dps
axh2 = mplot_dtor1r2d_addtitle(dps, axh2, opt);

return
m = mplot_signal_and_fit(S, xps, @dtor1r2d_1d_data2fit, @dtor1r2d_1d_fit2data, axh, opt);
title(axh,['\omega_{rms} / 2\pi = ' num2str(min(xps.momega(xps.momega>0)/2/pi),3) ' - ' num2str(max(xps.momega/2/pi),3) ' Hz'])

% Get dps from dtor1r2d_4d_fit2param
mfs.m = zeros(1,1,1,numel(m)); mfs.m(1,1,1,:) = m;
dps = dtor1r2d_4d_fit2param(mfs.m,[],opt);

[dpar,dperp,theta,phi,w] = dtor1r2d_4d_m2parso(mfs.m,mean(xps.momega));
[dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(squeeze(dpar),squeeze(dperp),squeeze(theta),squeeze(phi));
dtd_1x6 = [dxx dyy dzz sqrt(2)*[dxy dxz dyz]];

% Plot the tensor distribution
% [dtd_1x6,w] = dtd_dist2nx6w(dtd_m2dtd(m));
mplot_dtd(dtd_1x6, w, opt.dtor1r2d.dmin, opt.dtor1r2d.dmax, axh2, opt);
mplot_dtd_addstats(dps, axh2, opt);
mplot_dtd_addtitle(dps, axh2, opt);

