function m = fexi11_plot(S, xps, h, h2)

opt = fexi11_opt;

m = fexi11_1d_data2fit(S, xps, opt);

S_fit = fexi11_1d_fit2data(m, xps);


c_list = unique(xps.s_ind)';
tm      = zeros(1, numel(c_list));
ADC     = tm;
ADC_fit = tm;
bf      = tm;

for c = c_list
    
    ind = xps.s_ind == c;
    
    b = xps.b(ind);
    tm(c == c_list) = mean(xps.mde_tm12(ind));
    
    bf(c == c_list) = mean(xps.mde_b1(ind));
    
    p = polyfit(b, log(S(ind)), 1);
    ADC(c == c_list) = -p(1);
    
    p = polyfit(b, log(S_fit(ind)), 1);
    ADC_fit(c == c_list) = -p(1);
    
end

tm(bf <= 0.1e9) = -10e-3;

[~,ind] = sort(tm);
tm = tm(ind);
ADC = ADC(ind);
ADC_fit = ADC_fit(ind);

cla(h);
hold(h, 'off');
plot(h, tm * 1e3, ADC * 1e9, 'o', 'markerfacecolor', 'blue');
hold(h, 'on');
plot(h, tm * 1e3, ADC_fit * 1e9, 'k-');

ylabel(h, 'ADC (um^2/ms)');
xlabel(h, 'Mixing time (ms)');

xlim(h, [-20e-3 max(tm) + 20e-3] * 1e3);

title(h, sprintf('ADC = %0.1f, sigma = %0.2f, AXR = %0.1f', ...
    m(1) * 1e9, m(2), m(3)));

axis(h, 'on');
box(h, 'off');
set(h, 'tickdir','out');

cla(h2);
x = 1:numel(S);
plot(h2, x, S, 'bx', x, S_fit, 'k.-'); 

axis(h2, 'on');
box(h2, 'off');
set(h2, 'tickdir','out');

ylabel(h2, 'Signal^ _ ');
xlabel(h2, 'Volume^ _ ');


