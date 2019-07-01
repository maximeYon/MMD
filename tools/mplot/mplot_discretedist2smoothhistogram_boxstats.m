function axh = mplot_discretedist2smoothhistogram_boxstats(dist_d, axh, opt)
% function axh = mplot_discretedist2smoothhistogram_boxstats(dist_d, axh, opt)
%

if (nargin < 2), axhh = gca; end
if (nargin < 3), opt = []; end

opt = mplot_opt();

dist_s = dist_1d_discrete2boxplotpars(dist_d);
wmax = max(dist_s.w);

plot(axh,dist_s.x, dist_s.w,'b-','LineWidth',opt.mplot.lw)
hold(axh,'on')
plot(axh,dist_s.x,dist_s.percentile*max(dist_s.w),'k-','LineWidth',.5*opt.mplot.lw)

axis(axh,'tight')
set(axh,'Box','off','TickDir','out','TickLength',.02*[1 1],'YTick',[],'YLim',wmax*[-.1 1.2],'FontSize',opt.mplot.fs,'LineWidth',opt.mplot.lw)
title(axh,{['Q1 = ' num2str(dist_s.xQ1,3) ' mean = ' num2str(dist_s.xQ2,3) ' Q3 = ' num2str(dist_s.xQ3,3)];
    ['whisker down = ' num2str(dist_s.xwhiskerdown,3) ' up = ' num2str(dist_s.xwhiskerup,3)]}, 'FontWeight','normal')
boxheight = .1*wmax;
boxymid = 1.1*wmax;
boxx = [dist_s.xQ1*[1; 1]; dist_s.xQ3*[1; 1]; dist_s.xQ1];
boxy = boxymid + .5*boxheight*[-1 1 1 -1 -1]';
meanx = dist_s.xQ2*[1; 1];
meany = boxymid + .5*boxheight*[-1; 1];
whiskerx = [dist_s.xwhiskerdown*[1; 1] [dist_s.xwhiskerdown; dist_s.xQ1] [dist_s.xQ3; dist_s.xwhiskerup] dist_s.xwhiskerup*[1; 1]];
whiskery = boxymid + .5*boxheight*[[-1; 1] [0; 0] [0; 0] [-1; 1]];

plot(axh,meanx,meany,'k-','LineWidth',opt.mplot.lw)
plot(axh,boxx,boxy,'k-','LineWidth',opt.mplot.lw)
plot(axh,whiskerx,whiskery,'k-','LineWidth',opt.mplot.lw)

