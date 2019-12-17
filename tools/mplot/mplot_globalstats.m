function mplot_globalstats(method, dps, fig_fn, clim)
% function mplot_globalstats(method, dps, fig_fn, clim)

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

papersize = 17.5*[1 1/3];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

if ~isempty(fig_fn)
    msf_mkdir(fileparts(fig_fn));
    print(fig_fn,'-loose','-dpdf')
end



