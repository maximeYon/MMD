function m = gwf_maxwell_index(gwf, rf, dt)
% function m = gwf_maxwell_index(gwf, rf, dt)

wf_nx6 = tm_1x3_to_1x6(1, 0, gwf);
rf_nx6 = repmat(rf, 1, 6);

M = tm_1x6_to_3x3( sum( wf_nx6 .* rf_nx6 ) * dt );

m = sqrt(trace(M*M));

