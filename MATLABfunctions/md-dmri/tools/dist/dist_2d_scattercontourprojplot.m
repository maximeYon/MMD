function [hax,hscatter,hcontour,hprojx,hprojy] = dist_2d_scattercontourprojplot(hax,hax_projx,hax_projy,dist_d,contourpars,axpars)

lw_contour = 1;
lw_axes = 2;
lw_plot = 1;
fs_axes = 15;

axes(hax)
hold on

x = dist_d.x(:); x(~isfinite(x)) = 0;
y = dist_d.y(:); y(~isfinite(y)) = 0;
a = dist_d.a(:)+eps; a(~isfinite(a)) = 0;
r = dist_d.r(:);
g = dist_d.g(:);
b = dist_d.b(:);
c = repmat(dist_d.bright(:),[1 3]).*[r g b];

dist_d.n = numel(x);
dist_d.x = x;
dist_d.y = y;
dist_d.w = a;

dist_s.x = linspace(axpars.xmin,axpars.xmax,contourpars.Nx)';
dist_s.y = linspace(axpars.ymin,axpars.ymax,contourpars.Ny)';
dist_s.xsigma = 3*abs(dist_s.x(2) - dist_s.x(1));
dist_s.ysigma = 3*abs(dist_s.y(2) - dist_s.y(1));

dist_s = dist_2d_discrete2smooth(dist_d,dist_s);

C = contourc(double(dist_s.x),double(dist_s.y),double(dist_s.w'),contourpars.Nlevels);
hcontour = [];
count = 1;
while count < length(C)
    numxy = C(2,count);
    xtemp = C(1,count+(1:numxy));
    ytemp = C(2,count+(1:numxy));
    h = plot(xtemp,ytemp,'k-','LineWidth',lw_contour);
    hcontour = [hcontour; h];
    count = count + numxy + 1;
end

axis([axpars.xmin axpars.xmax axpars.ymin axpars.ymax])
set(hcontour,'Color',.5*[1 1 1])
set(hax,'LineWidth',lw_axes,'FontSize',fs_axes,'YAxisLocation','right','TickDir','out','Box','on')

wprojx = sum(double(dist_s.w),2);
wprojy = sum(double(dist_s.w),1);

%%% Mod Maxime 2022: projection normalization
Norm_fact = sum(wprojx);
wprojx = wprojx./Norm_fact;
wprojy = wprojy./Norm_fact;
wprojx(isnan(wprojx)) =0.0000001;wprojy(isnan(wprojy)) =0.0000001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hprojx = plot(hax_projx,dist_s.x,wprojx,'-');
hprojy = plot(hax_projy,wprojy,dist_s.y,'-');
set(hax_projx,'XLim',[axpars.xmin axpars.xmax],'YLim',max(abs(wprojx))*[-.1 1.5])
set(hax_projy,'YLim',[axpars.ymin axpars.ymax],'XLim',max(abs(wprojy))*[-.1 1.5],'XDir','reverse')
axis([hax_projx; hax_projy],'off')

if ~isfield(axpars,'no_scatter')
    hscatter = scatter(hax,x,y,a,c);
    set(hscatter,'LineWidth',lw_plot)
else
    hscatter = scatter(hax,0,0,eps,[1 1 1]);
end


