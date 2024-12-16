function axh = mplot_dtd_addstats(dps, axh, opt)
% function axh = mplot_dtd_addstats(dps, axh, opt)
%

hold(axh,'on')
x = log10(dps.mdiso);
y = log10((1+2*sqrt(dps.nmsdaniso))/(1-sqrt(dps.nmsdaniso)));
if ~isreal(y)
    y = abs(y);
    if any([y > log10(opt.dtd.dmax/opt.dtd.dmin) dps.nmsdaniso > 1]);
        y = log10(opt.dtd.dmax/opt.dtd.dmin);
    end
end
if isfield(dps,'signaniso')
    if dps.signaniso < 0
        y = log10((1-2*sqrt(dps.nmsdaniso))/(1+sqrt(dps.nmsdaniso)));
    end
    if ~isreal(y)
        y = -log10(opt.dtd.dmax/opt.dtd.dmin);
    end
end
x_low = log10(dps.mdiso-sqrt(dps.vdiso));
x_high = log10(dps.mdiso+sqrt(dps.vdiso));
if ~isreal(x_low)
    x_low = abs(x_low);
    x_high = abs(x_high);
end

ph = plot(axh,x,y,'ko',[x_low x_high],y*[1 1],'k-');
set(ph,'Color',.7*[1 1 1],'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',.7*[1 1 1])
