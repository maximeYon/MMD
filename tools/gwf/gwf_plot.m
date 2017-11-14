function gwf_plot(gwf, rf, dt, col)
% function gwf_plot(gwf, rf, dt, col)
%
% gwf - n x (1, 2, or 3)  [gradient as executed by amplifier]
% rf  - n x 1             [effect of RF pulses, either 1 or -1]
% dt  - 1 x 1             [time step]
% col - 1 x 3 cell        [colours]

if (nargin < 2), error('rf needed'); end
if (nargin < 3), error('dt needed'); end
if (nargin < 4), col = {}; end

% init and test things
if (isempty(rf)), rf = gwf_to_rf(gwf); end
if (numel(col) < 3), col = {[0.2 0.5 0.95], [0.4 0.9 0.1], [0.99 0.3 0.2]}; end

gwf_check(gwf, rf, dt);    

if (size(gwf, 2) == 1)
    gwf = cat(2, gwf, zeros(size(gwf)), zeros(size(gwf)));
end

% compute some things
g_eff = gwf .* repmat(rf, 1, 3);
t = linspace(0, size(gwf, 1) * dt, size(gwf, 1));

for c = 1:size(gwf, 2)
    area(t * 1e3, g_eff(:, c) * 1e3, ...
        'FaceAlpha', 0.7, ...
        'EdgeColor', 'black', ...
        'FaceColor', col{c}, ...
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






