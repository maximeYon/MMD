function mplot_globalstats(method, dps, fig_fn, clim, opt)
% function mplot_globalstats(method, dps, fig_fn, clim, opt)

%Plot parameter maps
figure(2), clf

Nparam = 3;
left = .05;
bottom = .17;
dleft = (1-left)/(Nparam+0);
width = .9*dleft;
height = .7;

opt = mplot_opt();
opt.mplot.fs = 7;

dist_d.w = dps.s0(:);
dist_d.x = dps.mdiso(:);
axh = axes('position',[left bottom width height]);
axh = mplot_discretedist2smoothhistogram_boxstats(dist_d, axh, opt);
xlabel('size, E[D_{iso}] / m^2s^{-1}')
xlim = clim.mdiso + .05*abs(diff(clim.mdiso))*[-1 1];
set(gca,'XLim',xlim)

left = left+dleft;
dist_d.x = dps.msddelta(:);
axh = axes('position',[left bottom width height]);
axh = mplot_discretedist2smoothhistogram_boxstats(dist_d, axh, opt);
xlabel('shape, E[D_\Delta^2]')
xlim = clim.msddelta + .05*abs(diff(clim.msddelta))*[-1 1];
set(gca,'XLim',xlim)

if strcmp(method,'dtr2d')
    left = left+dleft;
    dist_d.x = dps.mr2(:);
    axh = axes('position',[left bottom width height]);
    axh = mplot_discretedist2smoothhistogram_boxstats(dist_d, axh, opt);
    xlabel('E[R_2] / s^{-1}')
    xlim = clim.mr2 + .05*abs(diff(clim.mr2))*[-1 1];
    set(gca,'XLim',xlim)
elseif strcmp(method,'dtr1d')
    left = left+dleft;
    dist_d.x = dps.mr1(:);
    axh = axes('position',[left bottom width height]);
    axh = mplot_discretedist2smoothhistogram_boxstats(dist_d, axh, opt);
    xlabel('E[R_1] / s^{-1}')
    xlim = clim.mr1 + .05*abs(diff(clim.mr1))*[-1 1];
    set(gca,'XLim',xlim)
end

if strcmp(method,'dtd')
    opt = mplot_opt(opt);

    left = left+dleft;
    dist_d.w = dps.s0(:);
    dist_d.x = dps.mdiso(:);
    dist_d.y = dps.msddelta(:);
    dist_d.n = numel(dist_d.w);

    axpars.xmin = clim.mdiso(1);
    axpars.xmax = clim.mdiso(2);
    axpars.ymin = clim.msddelta(1);
    axpars.ymax = clim.msddelta(2);
    contourpars.Nx = 50;
    contourpars.Ny = contourpars.Nx;
    contourpars.Nlevels = 20;

    dist_s.x = linspace(axpars.xmin,axpars.xmax,contourpars.Nx)';
    dist_s.y = linspace(axpars.ymin,axpars.ymax,contourpars.Ny)';
    dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));
    dist_s.ysigma = 2*abs(dist_s.y(2) - dist_s.y(1));

    dist_s = dist_2d_discrete2smooth(dist_d,dist_s);

    C = contourc(double(dist_s.x),double(dist_s.y),double(dist_s.w'),contourpars.Nlevels);

    axh = axes('position',[left bottom width height]);
    hold(axh,'on')

    hcontour = [];
    count = 1;
    while count < length(C)
        numxy = C(2,count);
        xtemp = C(1,count+(1:numxy));
        ytemp = C(2,count+(1:numxy));
        h = plot(axh,xtemp,ytemp,'k-','LineWidth',.5*opt.mplot.lw);
        hcontour = [hcontour; h];
        count = count + numxy + 1;
    end

    xlabel('size, E[D_{iso}] / m^2s^{-1}')
    xlim = clim.mdiso + .05*abs(diff(clim.mdiso))*[-1 1];
    ylabel('shape, E[D_\Delta^2]')
    ylim = clim.msddelta + .05*abs(diff(clim.msddelta))*[-1 1];
    set(gca,'XLim',xlim)
    set(gca,'YLim',ylim)
    axis square
    set(axh,'Box','off','TickDir','out','TickLength',.02*[1 1],'FontSize',opt.mplot.fs,'LineWidth',opt.mplot.lw)

end

papersize = 17.5*[1 1/3];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

if ~isempty(fig_fn)
    msf_mkdir(fileparts(fig_fn));
    print(fig_fn,'-loose','-dpdf')
end



