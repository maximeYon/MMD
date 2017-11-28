function [gwf_out, rf_out, dt_out] = gwf_interpolate(gwf, rf, dt, scale)
% function [gwf, rf, dt] = gwf_interpolate(gwf, rf, dt, scale)
%
% Interpolate the gwf -- currently nearest neighbour
%
% scale -- integer

gwf_check(gwf, rf, dt);

scale = round(scale);

gwf_out = zeros(size(gwf,1) * scale, 3);
rf_out  = zeros(size(rf, 1) * scale, 1);

for k = 1:scale
    gwf_out( ((1:size(gwf,1))-1) * scale + k, :) = gwf;
    rf_out ( ((1:size(gwf,1))-1) * scale + k, :) = rf;
end

dt_out = dt / scale;

