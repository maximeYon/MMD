function dtd_plot(S, xps, axh, axh2)
% function dtd_plot(S, xps, axh, axh2)

if (nargin < 4), axh2 = []; end

ms = 5;
fs = 10;
lw = 1;

opt = mdm_opt();
opt = dtd_opt(opt);

opt.dtd.dmin = .2/max(xps.b);

m = feval('dtd_1d_data2fit', S, xps, opt);
S_fit = feval('dtd_1d_fit2data', m, xps)';

cla(axh);
hold(axh, 'off');
h1 = plot(axh,1:xps.n,S,'o',1:xps.n,S_fit,'x');
set(h1,'MarkerSize',ms,'LineWidth',lw)
axis(axh,'tight')
set(axh,'XLim',xps.n*[-.1 1.1], 'YLim',max(S_fit)*[-.1 1.1],...
    'Box','off','TickDir','out','TickLength',.02*[1 1],...
'FontSize',fs,'LineWidth',lw)
xlabel(axh,'Acq number','FontSize',fs)
ylabel(axh,'Signal','FontSize',fs)


dmin = opt.dtd.dmin;
dmax = opt.dtd.dmax;
ratiomin = dmin/dmax;
ratiomax = dmax/dmin;

xmin = log10(dmin);
xmax = log10(dmax);
ymin = log10(1/ratiomax);
ymax = log10(ratiomax);

[n,par,perp,theta,phi,w] = dtd_dist2par(dtd_m2dtd(m));
s0 = sum(w);

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
    c.ms = 5*ms*sqrt(w/s0);
    c.bright = fa;
    c.r = abs(xcos);
    c.g = abs(ycos);
    c.b = abs(zcos);

    for nc = 1:n
        h1 = plot(axh2,c.x(nc),c.y(nc),'o','LineWidth',.01);
        hold on
        col = c.bright(nc)*[c.r(nc) c.g(nc) c.b(nc)];
        set(h1,'MarkerSize',c.ms(nc),'Color',col,'MarkerFaceColor',col)
    end
end

set(axh2,'XLim',[xmin xmax], 'YLim',[ymin ymax],'YAxisLocation','right',...
'XTick',[-11:.5:-8],'YTick',-2:.5:2,'TickDir','out','TickLength',.02*[1 1],...
'FontSize',fs,'LineWidth',lw,'Box','on')
axis(axh2,'square')
xlabel(axh2,'size, log(\itD\rm_{iso} / m^2s^-^1)','FontSize',fs)
ylabel(axh2,'shape, log(\itD\rm_{||} / \itD\rm_{\perp})','FontSize',fs)
title(axh2,['orientation, [RGB]=[xyz]'],'FontSize',fs,'FontWeight','normal')
