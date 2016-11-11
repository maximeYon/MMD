function dtd_pa_plot(m, opt)
% function dtd_pa_mkpdf(dps_fn, pdf_path, opt)


dmin = opt.dtd_pa.dmin;
dmax = opt.dtd_pa.dmax;
ratiomax = dmax/dmin;

nx = 50; ny = nx;

sigma_x = .2;
sigma_y = sigma_x;

xmin = log10(dmin)-3*sigma_x;
xmax = log10(dmax)+3*sigma_x;
ymin = log10(1/ratiomax)-1*sigma_y;
ymax = log10(ratiomax)+1*sigma_y;

x = linspace(xmin,xmax,nx)';
dx = x(2)-x(1);
y = linspace(ymin,ymax,ny)';
dy = y(2)-y(1);

[xx,yy] = ndgrid(x,y);

np = nx*ny;

ptot = [];
for c = 1:size(m,1)
    dtd = dtd_pa_m2dtd(m(c,:));
    [n,par,perp,w] = dtd_pa_dist2par(dtd);
    
    if (n == 0), continue; end;
    
    logiso = log10((par + 2*perp)/3);
    logratio = log10(par./perp);
    
    x0 = logiso(:);
    y0 = logratio(:);
    nx0 = n;
    xx_k = repmat(xx(:),[1 nx0]);
    yy_k = repmat(yy(:),[1 nx0]);
    x0_k = repmat(x0',[nx*ny 1]);
    y0_k = repmat(y0',[nx*ny 1]);
    
    k = 1/(sigma_x*sqrt(2*pi)).*exp(-(xx_k-x0_k).^2/(2*sigma_x^2)).*...
        1/(sigma_y*sqrt(2*pi)).*exp(-(yy_k-y0_k).^2/(2*sigma_y^2));
    
    if (isempty(ptot))
        ptot = k*w;
    else
        ptot = ptot + k*w;
    end
end

z = reshape(ptot,[nx ny]);

zprojx = sum(z,2);
zprojy = sum(z,1);

nclevels = 3;
clevels = max(z(:))*linspace(0,1,nclevels+2);
clevels = clevels(2:(nclevels+1));

left    = .6;
bottom  = .5;
height  = .4;
width   = 0.3; 
proj_height = height * 0.2;
proj_width = width * 0.2;
lw = 1;
fs = 12;

clf;

axes('position',[left bottom width height])
contour(x,y,z',clevels,'k','LineWidth',.5*lw)
set(gca,'XLim',[xmin xmax], 'YLim',[ymin ymax],'YAxisLocation','right',...
    'XTick',[-11:-8],'YTick',-2:2,'TickDir','out','TickLength',.03*[1 1],...
    'FontSize',fs,'LineWidth',lw)
% axis square

xlabel('log(D_{iso} / m^2s^-^1)','FontSize',fs)
ylabel('log(D_{||} / D_{\perp})','FontSize',fs)

axes('position',[left bottom+height width proj_height])
plot(x,zprojx/max(zprojx),'k-')
set(gca,'XLim',[xmin xmax], 'YLim',[-.3 1.1])
axis off

axes('position',[left-proj_width bottom proj_width height])
plot(zprojy/max(zprojy),y,'k-')
set(gca,'YLim',[ymin ymax], 'XLim',[-.3 1.1],'XDir','reverse')
axis off







