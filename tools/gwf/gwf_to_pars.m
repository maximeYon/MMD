function xps = gwf_to_pars(gwf, rf, dt)
% function xps = gwf_to_pars(gwf, rf, dt)
%
% gwf - gradient waveform of size N x 3
% rf  - effect of rf pulses (range -1 to 1), size N x 1
% dt  - time step of waveform

if (isempty(rf)), rf = gwf_to_rf(gwf); end

if (iscell(gwf))
    xps = cell(numel(gwf));
    for c = 1:numel(gwf)
        xps{c} = gwf_to_pars(gwf{c}, rf, dt);
    end
    xps = mdm_xps_merge(xps);
    xps = rmfield(xps, 's_ind');
    return;
end


% Effective gradient waveform,  time, and other useful functions
t = (1:size(gwf,1))' * dt;
r = @(rf,n) repmat(rf, 1, n);
g_eff = gwf .* r(rf, 3);



% b-tensor and associated metrics
xps = mdm_xps_from_bt(gwf_to_bt(gwf, rf, dt));

% Max gradient amplitude and slew rate 
xps.slew_rate_max     = max(abs(diff(gwf, 1, 1)),[],1) / dt;
xps.gradient_max      = max(abs(gwf),[],1);
xps.gradient_max_norm = max(abs( sqrt(sum(gwf.^2, 2))));

% Spectral properties
[M0,M1] = gwf_to_spectral_moments(gwf, rf, dt);
M1_pars = tm_3x3_to_tpars(M1);
xps.gwf_spectral_M0 = tm_3x3_to_1x6(M0);
xps.gwf_spectral_M1 = tm_3x3_to_1x6(M1);
xps.gwf_spectral_M1_iso   = M1_pars.iso;
xps.gwf_spectral_M1_delta = M1_pars.delta;

% Approximately V_omega^1/2 in Nilsson et al (2017) NMR Biomed
xps.gwf_spectral_trace_rms_norm = ...
    sqrt(xps.gwf_spectral_M1_iso / xps.b);



% Maxwell index: This is a relative index
M = tm_1x6_to_3x3( sum(tm_1x3_to_1x6(1, 0, gwf) .* r(rf, 6)) * dt );

xps.maxwell_index = 1e9 * sqrt(trace(M * M));



% Energy index
xps.energy_index = sum(gwf(:).^2) * dt * 1e4;
 


% Gradient moments
xps.k0 = msf_const_gamma('1H') * sum( g_eff .* r(t.^0,3), 1 ) * dt;
xps.k1 = msf_const_gamma('1H') * sum( g_eff .* r(t.^1,3), 1 ) * dt;
xps.k2 = msf_const_gamma('1H') * sum( g_eff .* r(t.^2,3), 1 ) * dt;

xps.alpha_norm = (sum(xps.k1.^2));



