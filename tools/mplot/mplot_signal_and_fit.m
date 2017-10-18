function m = mplot_signal_and_fit(S, xps, f_data2fit, f_fit2data, h, opt)
% function mplot_signal_and_fit(S, xps, f_data2fit, f_fit2data, opt)

opt = mplot_opt(opt);

% Fit and predict signal
m     = f_data2fit(S, xps, opt); 
S_fit = f_fit2data(m, xps)'; 

% Clear, plot, and configure 
cla(h); hold(h, 'off');

plot(h,1:xps.n,S,    'ko','MarkerSize',round(opt.mplot.ms*1.5), 'MarkerFaceColor', 'black');
hold(h, 'on');
plot(h,1:xps.n,S_fit,'ko','MarkerSize',round(opt.mplot.ms*1.0), 'MarkerFaceColor', 'red');
hold(h, 'off');

axis(h,'tight');

set(h,...
    'XLim',xps.n*[-.1 1.1], ...
    'YLim',max(S_fit)*[-.1 1.1],...
    'Box','off',...
    'TickDir','out',...
    'TickLength',.02*[1 1],...
    'FontSize',opt.mplot.fs,...
    'LineWidth',opt.mplot.lw);

xlabel(h,'Acq number','FontSize',opt.mplot.fs)
ylabel(h,'Signal','FontSize',opt.mplot.fs)