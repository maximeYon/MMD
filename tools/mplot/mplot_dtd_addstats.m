function axh = mplot_dtd_addstats(dps, axh, opt)
% function axh = mplot_dtd_addstats(dps, axh, opt)
%

hold(axh,'on')
x = log10(dps.mdiso);
y = log10((1+2*sqrt(dps.msdanison))/(1-sqrt(dps.msdanison)));
if ~isreal(y)
    y = log10(opt.dtd.dmax/opt.dtd.dmin);
end
if isfield(dps,'signaniso')
    if dps.signaniso < 0
        y = log10((1-2*sqrt(dps.msdanison))/(1+sqrt(dps.msdanison)));
    end
    if ~isreal(y)
        y = -log10(opt.dtd.dmax/opt.dtd.dmin);
    end
end
x_low = log10(dps.mdiso-sqrt(dps.vdiso));
x_high = log10(dps.mdiso+sqrt(dps.vdiso));
ph = plot(axh,x,y,'ko',[x_low x_high],y*[1 1],'k-');
set(ph,'Color',.7*[1 1 1],'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',.7*[1 1 1])
