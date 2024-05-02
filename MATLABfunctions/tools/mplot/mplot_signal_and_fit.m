function m = mplot_signal_and_fit(S, xps, f_data2fit, f_fit2data, h, opt)
% function mplot_signal_and_fit(S, xps, f_data2fit, f_fit2data, opt)

opt = mplot_opt(opt);

% Fit and predict signal
m     = f_data2fit(S, xps, opt); 
S_fit = f_fit2data(m, xps)'; 

% Clear, plot, and configure 
cla(h); hold(h, 'off');

if isfield(xps,'theta')
    xps_array = [round(xps.b,-7), round(xps.b_delta,1), xps.theta/pi*180, xps.phi/pi*180];
    [xps_array_sort, ind_sort] = sortrows(xps_array,[1 2 3 4]);
else
    xps_array = [round(xps.b,-7), round(xps.b_delta,1)];
    [xps_array_sort, ind_sort] = sortrows(xps_array,[1 2]);
end

plot(h,1:xps.n,S(ind_sort),    'ro','MarkerSize',round(opt.mplot.ms*.5), 'MarkerFaceColor', 'none');
hold(h, 'on');
%plot(h,1:xps.n,S_fit,'ro','MarkerSize',round(opt.mplot.ms*1.0), 'MarkerFaceColor', 'red');
plot(h,1:xps.n,S_fit(ind_sort),'k.','MarkerSize',round(opt.mplot.ms*1.0));
hold(h, 'off');

axis(h,'tight');

set(h,...
    'XLim',xps.n*[-.1 1.1], ...
    'YLim',max([S_fit(:); S(:)])*[-.1 1.1],...
    'Box','off',...
    'TickDir','out',...
    'TickLength',.02*[1 1],...
    'FontSize',opt.mplot.fs,...
    'LineWidth',opt.mplot.lw);

xlabel(h,'Acq number','FontSize',opt.mplot.fs)
ylabel(h,'Signal','FontSize',opt.mplot.fs)