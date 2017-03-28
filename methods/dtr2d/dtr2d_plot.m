function dtr2d_plot(S, xps, axh, axh2)
% function dtr2d_plot(S, xps, axh, axh2)

if (nargin < 4), axh2 = []; end

ms = 5;
fs = 10;
lw = 1;

opt = mdm_opt();
opt = dtr2d_opt(opt);

opt.dtr2d.dmin = .1/max(xps.b);
opt.dtr2d.r2min = .1/max(xps.te);
opt.dtr2d.r2max = 1/min(xps.te);

m = feval('dtr2d_1d_data2fit', S, xps, opt);
S_fit = feval('dtr2d_1d_fit2data', m, xps)';


dmin = opt.dtr2d.dmin;
dmax = opt.dtr2d.dmax;
ratiomin = dmin/dmax;
ratiomax = dmax/dmin;
r2min = opt.dtr2d.r2min;
r2max = opt.dtr2d.r2max;

xmin = log10(dmin);
xmax = log10(dmax);
ymin = log10(1/ratiomax);
ymax = log10(ratiomax);
zmin = log10(r2min);
zmax = log10(r2max);

[n,par,perp,theta,phi,r2,w] = dtr2d_dist2par(dtr2d_m2dtr2d(m));
s0 = sum(w);

cla(axh);
hold(axh, 'off');
h1 = plot(axh,1:xps.n,S,'o',1:xps.n,S_fit,'x');
set(h1,'MarkerSize',ms,'LineWidth',lw)
axis(axh,'tight')
set(axh,'XLim',xps.n*[-.1 1.1], 'YLim',s0*[-.1 1.1],...
    'Box','off','TickDir','out','TickLength',.02*[1 1],...
'FontSize',fs,'LineWidth',lw)
xlabel(axh,'Acq number','FontSize',fs)
ylabel(axh,'Signal','FontSize',fs)

cla(axh2);
hold(axh2, 'on');

if n>0
    xcos = cos(phi).*sin(theta);
    ycos = sin(phi).*sin(theta);
    zcos = cos(theta);

    iso = tm_eigvals2iso([par perp perp]);
    fa = tm_eigvals2fa([par perp perp]);

    c.x = log10(iso);
    c.y = log10(par./perp);
    c.z = log10(r2) - log10(r2min)+log10(ratiomax);
    c.ms = 5*ms*sqrt(w/s0);
    c.bright = fa;
    c.r = abs(xcos);
    c.g = abs(ycos);
    c.b = abs(zcos);

    for nc = 1:n
        h1 = plot(axh2,c.x(nc),c.y(nc),'o','LineWidth',.01);
        h2 = plot(axh2,c.x(nc),c.z(nc),'o','LineWidth',.01);
%        hold on
        col = c.bright(nc)*[c.r(nc) c.g(nc) c.b(nc)];
        set([h1; h2],'MarkerSize',c.ms(nc),'Color',col,'MarkerFaceColor',col)
    end
end
plot(axh2,[log10(dmin) log10(dmax)],log10(ratiomax)*[1 1],'k-','LineWidth',lw);

set(axh2,'XLim',[xmin xmax], 'YLim',[ymin ymax+zmax-zmin],'YAxisLocation','right',...
'XTick',[-11:.5:-8],'YTick',-2:.5:2,'TickDir','out','TickLength',.02*[1 1],...
'FontSize',fs,'LineWidth',lw,'Box','on')
axis(axh2,'square')
xlabel(axh2,'size, log(\itD\rm_{iso} / m^2s^-^1)','FontSize',fs)
ylabel(axh2,'shape, log(\itD\rm_{||} / \itD\rm_{\perp})    log(\it{R}\rm_2 / s^-^1)','FontSize',fs)
title(axh2,['orientation, [RGB]=[xyz]'],'FontSize',fs,'FontWeight','normal')
