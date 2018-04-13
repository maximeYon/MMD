function m = gwf_maxwell_index(gwf, rf, dt)
% function m = gwf_maxwell_index(gwf, rf, dt)
%
% Compute the Maxwell index from a gradient waveform

if (isempty(rf)), rf = gwf_to_rf(gwf); end

wf_Nx6 = tm_1x3_to_1x6(1, 0, gwf);
rf_Nx6 = repmat(rf, 1, 6);

M = tm_1x6_to_3x3( sum( wf_Nx6 .* rf_Nx6 ) * dt );

m = sqrt(trace(M*M));

