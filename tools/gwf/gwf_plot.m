function gwf_plot(gwf, dt)
% function gwf_plot(gwf, dt)
%
% gwf - (1, 2, or 3) x n

if (size(gwf,1) > 3), error('gwf should be (1, 2, or 3) x n'); end

for c = 1:size(gwf, 1)
    if (abs(sum(gwf(c,:))) > 1e-3)
        warning('gwf does not integrate to zero');
    end
end

t = linspace(0, size(gwf, 2) * dt, size(gwf, 2));

col = {[0.9 0.1 0.0], [0.3 0.8 0.0], [0.1 0.4 0.9]};

clf;
set(gcf,'color','white');

% PLOT THE GRADIENT WAVEFORM
subplot(2,2,1);
for c = 1:size(gwf, 1)
    plot(t * 1e3, gwf(c, :) * 1e3, 'color', col{c}, 'linewidth', 2); hold on;
end

xlabel('t [ms]');
ylabel('g [mT/m]');

xlim([min(t) max(t)] * 1e3);
ylim([-1 1] * max(abs(gwf(:) * 1e3)));

box off;
set(gca,'tickdir', 'out', 'ticklength', [0.03 0.05]);
set(gca,'linewidth', 2);

% PLOT THE ASSOCIATED SPECTRA
subplot(2,2,3);
vf = zeros(size(gwf,1), 1);
b = vf;
for c = 1:size(gwf, 1)
    
    % zero fill before power computation
    tmp = gwf(c,:);
    tmp = cat(2, zeros(size(tmp)), tmp, zeros(size(tmp)));
    tmp = cat(2, zeros(size(tmp)), tmp, zeros(size(tmp)));
    tmp = cat(2, zeros(size(tmp)), tmp, zeros(size(tmp)));
    
    [ps, f, df] = gwf_power_spectrum_from_g(tmp, dt);
        
    
    % compute the second moment
    b(c) = gwf_b_from_g(tmp, dt);
    vf(c) = gwf_vf_from_g(tmp, dt);
    plot(f, ps * 1e-9, 'color', col{c}, 'linewidth', 2); hold on;
end

xlim([-1 1] * 5 * sqrt(max(vf)));
xlabel('f [Hz]');
ylabel('Enc power'); 
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.03 0.05]);
set(gca,'linewidth', 2);
set(gca,'ytick', []);


% COMPUTE THINGS FOR ANALYSIS
slew_rate = max(abs(diff(gwf, 1, 2)),[],2) / dt; 

subplot(1,2,2);
text(0,0,{...
    'GRADIENT WAVEFORM ANALYSIS', ...
    ['max slew rate = ' num2str(slew_rate',3) ' T/m/s'], ...
    ['b = ' num2str(b' * 1e-9,2) ' ms/um^2'], ...
    ['V_f = ' num2str(vf', 3) ' 1/s^2'], ...
    ['"Diffusion time" = ' num2str(1e3./vf',3) ' ms']});
ylim([-100 10]);
axis off;

