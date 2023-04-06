function [g_eff, t] = gwf_to_g_eff(gwf, rf, dt, opt)
% function gwf_to_g_eff(gwf, rf, dt, opt)
%
% Compute the effective gradient from the physical gradient


if (isempty(rf)), rf = gwf_to_rf(gwf); end
if (size(gwf, 2) == 1), gwf = cat(2, gwf, zeros(numel(gwf), 2)); end

gwf_check(gwf, rf, dt);

g_eff = gwf .* repmat(rf, 1, 3);
t = ((1:numel(rf))-1)' * dt; 


