function gwf_analysis(gwf, rf, dt, col)
% function gwf_analysis(gwf, rf, dt, col)
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
pars = gwf_to_pars(gwf, rf, dt);
[gwf_ps, f, df] = gwf_power_spectrum(gwf, rf, dt);
g_eff = gwf .* repmat(rf, 1, 3);

% set up axis
t = linspace(0, size(gwf, 1) * dt, size(gwf, 1));
p.f_max = 200;

clf;
set(gcf,'color','white');

% PLOT THE GRADIENT WAVEFORM
% --------------------------
subplot(2,2,1);
gwf_plot(gwf, rf, dt, col);

box off;
set(gca,'tickdir', 'out', 'ticklength', [0.03 0.05]);
set(gca,'linewidth', 2);
title('Effective gradient');

% PLOT THE SPECTRA
% --------------------------
subplot(2,2,3);
for c = 1:size(gwf, 2)
    
    area(f, gwf_ps(:,c) * 1e-9, ...
        'FaceAlpha', 0.7, ...
        'EdgeColor', 'black', ...
        'FaceColor', col{c}, ...
        'LineWidth', 2); hold on;
end

xlim([0 1] * max(p.f_max, 2 * pars.gwf_spectral_trace_rms_norm));
xlabel('f [Hz]');
ylabel('Enc power');
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.03 0.05]);
set(gca,'linewidth', 2);
set(gca,'ytick', []);
title('Power spectrum of q(t)');

% SHOW g-SPECTRA FOR MECHANICAL ANALYSIS
% --------------------------------------
subplot(2,2,4);

for c = 1:size(gwf, 2)
    
    % zero fill before power computation
    tmp = gwf(:, c);
    tmp = cat(1, zeros(size(tmp)), tmp, zeros(size(tmp)));
    tmp = cat(1, zeros(size(tmp)), tmp, zeros(size(tmp)));
    tmp = cat(1, zeros(size(tmp)), tmp, zeros(size(tmp)));
    
    ps = abs(fftshift(fft(tmp * dt))).^2;
    f  = linspace(-1/dt, 1/dt, numel(ps)) / 2;
    
    area(f, ps * 1e-9, ...
        'FaceAlpha', 0.7, ...
        'EdgeColor', 'black', ...
        'FaceColor', col{c}, ...
        'LineWidth', 2); hold on;
end

xlim([0 1] * max(p.f_max, 5 * pars.gwf_spectral_trace_rms_norm));
xlabel('f [Hz]');
ylabel('Enc power');
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.03 0.05]);
set(gca,'linewidth', 2);
set(gca,'ytick', []);
title('Power spectrum of g(t)');

% PRINT AN ANALYSIS
% --------------------------

axes('position', [0.54 0.5 0.4 0.35]);
text(0,0,{...
    'GRADIENT WAVEFORM ANALYSIS^ _ ', ...
    [sprintf('|k_0| = %2.2f 1/mm', 1e-3 * sqrt(sum(pars.k0.^2)))], ...
    ['b = ' sprintf('%2.2f', pars.b * 1e-9) ' ms/um^2, ' ...
    'b_\Delta = ' sprintf('%1.2f', pars.b_delta) ', ' ...
    'b_\eta = '   sprintf('%1.2f', pars.b_eta) '_ ^ '], ...    
    sprintf('|\\alpha|^2 = %1.0f s^2/mm^2 (flow comp)_ ^ ', pars.alpha_norm * 1e-6), ...
    sprintf('V_f = Tr(M_1) / Tr(b) =  %1.0f 1/s, [M_1]_\\Delta = %1.1f _ ^ ', pars.gwf_spectral_trace_rms_norm, pars.gwf_spectral_M1_delta), ...
    ['Max slew rate = ' sprintf('%2.0f, %2.0f, %2.0f', pars.slew_rate_max(1), pars.slew_rate_max(2), pars.slew_rate_max(3)) ' T/m/s^ _ '], ...
    ['Max gradient = ' sprintf('%2.0f, %2.0f, %2.0f', pars.gradient_max(1)*1e3, pars.gradient_max(2)*1e3, pars.gradient_max(3)*1e3) ' mT/m^ _ '], ...
    ['Energy index = ' sprintf('%2.0f', pars.energy_index) '_ ^ '], ...
    ['Maxwell index = ' sprintf('%2.0f', pars.maxwell_index * 1e9) ' (mT/m)^2 ms _ ^ '], ...
    });
ylim([-20 10]);
axis off;


%     ['V_f = ' num2str(vf', 4) ' 1/s^2_ '], ...
%     ['<V_f> = ' num2str(mean(vf'), 4) ' 1/s^2_ '], ...
%     ['"Diffusion time" = ' num2str(1e3./vf',3) ' ms^ _ '], ...







