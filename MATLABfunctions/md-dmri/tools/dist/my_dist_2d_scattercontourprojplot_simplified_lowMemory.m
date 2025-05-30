function [hax,hscatter,hcontour,hprojx,hprojy,wprojx,wprojy] = my_dist_2d_scattercontourprojplot_simplified_lowMemory(hax,hax_projx,hax_projy,dist_d,contourpars,axpars,max_projXY,my_yticks,xticks_yes,yticks_yes,my_sigma)

lw_contour = 1;
lw_axes = 1;
lw_plot = 1;
fs_axes = 8;

axes(hax)
hold on

x = dist_d.x(:); x(~isfinite(x)) = 0;
y = dist_d.y(:); y(~isfinite(y)) = 0;
a = dist_d.a(:)+eps; a(~isfinite(a)) = 0;

dist_d.n = numel(x);
dist_d.x = x;
dist_d.y = y;
dist_d.w = a;

dist_s.x = linspace(axpars.xmin,axpars.xmax,contourpars.Nx)';
dist_s.y = linspace(axpars.ymin,axpars.ymax,contourpars.Ny)';
dist_s.xsigma = my_sigma*abs(dist_s.x(2) - dist_s.x(1));
dist_s.ysigma = my_sigma*abs(dist_s.y(2) - dist_s.y(1));

dist_s = dist_2d_discrete2smooth_lowMemory(dist_d,dist_s);

C = contourc(double(dist_s.x),double(dist_s.y),double(dist_s.w'),contourpars.Nlevels);
axis square;
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
if strcmp(xticks_yes,'Yes');xticks([0 1 2 3 4 5]);else;xticks([]);end
% if strcmp(xticks_yes,'Yes');xticks([0 0.5 1 1.5 2 2.5 3 3.5 4 4.5]);else;xticks([]);end
if strcmp(yticks_yes,'Yes');yticks(my_yticks);else;yticks([]);end

wprojx = sum(double(dist_s.w),2);
wprojy = sum(double(dist_s.w),1);

%%% Mod Maxime 2022: projection normalization
Norm_fact = sum(wprojx);
wprojx = wprojx./Norm_fact;
wprojy = wprojy./Norm_fact;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hprojx = plot(hax_projx,dist_s.x,wprojx,'-');
hprojy = plot(hax_projy,wprojy,dist_s.y,'-');
if isempty(max_projXY)
    max_projX = max(abs(wprojx));
    max_projY = max(abs(wprojy)); 
else
    max_projX = max([max_projXY(1,1) max(abs(wprojx))]);
    max_projY = max([max_projXY(1,2) max(abs(wprojy))]);
end
set(hax_projx,'XLim',[axpars.xmin axpars.xmax],'YLim',max_projX*[-.1 1.8])
set(hax_projy,'YLim',[axpars.ymin axpars.ymax],'XLim',max_projY*[-.1 1.8],'XDir','reverse')
axis([hax_projx; hax_projy],'off')
set(hax_projx,'LineWidth',lw_axes,'FontSize',fs_axes)
set(hax_projy,'LineWidth',lw_axes,'FontSize',fs_axes)

if ~isfield(axpars,'no_scatter')
    hscatter = scatter(hax,x,y,a,c);
    set(hscatter,'LineWidth',lw_plot)
else
    hscatter = scatter(hax,0,0,eps,[1 1 1]);
end


