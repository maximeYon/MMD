function gwf_plot(gwf, rf, dt, opt)
% function gwf_plot(gwf, rf, dt, opt)
%
% gwf - n x (1, 2, or 3)  [gradient as executed by amplifier]
% rf  - n x 1             [effect of RF pulses, either 1 or -1]
% dt  - 1 x 1             [time step]
% opt - option structure

if (nargin < 4), opt = []; end

opt = gwf_opt(opt);

% init and test things
[g_eff, t] = gwf_to_g_eff(gwf, rf, dt);


for c = 1:size(gwf, 2)
    area(t * 1e3, g_eff(:, c) * 1e3, ...
        'FaceAlpha', 0.7, ...
        'EdgeColor', 'black', ...
        'FaceColor', opt.gwf.col{c}, ...
        'LineWidth', 2); hold on;
end

xlabel('t [ms]');
ylabel('g [mT/m]');

xlim([min(t - 3e-3) max(t + 3e-3)] * 1e3);
ylim([-1 1] * max(abs(gwf(:) * 1e3) * 1.1));

box off;
set(gca,'tickdir', 'out', 'ticklength', [0.03 0.05]);
set(gca,'linewidth', 2);
title('Effective gradient');






