function s = dnp_d1d8p15_1d_fit2data(m, xps, opt)
% function s = dnp_d1d8p15_1d_fit2data(m, xps, opt)
%
% m = [s0 psd rdnp rsd rcp r1 r1rho]

if nargin < 3
    opt = [];
end
    
% convert to readable parameters
s0 = m(1); % global signal max
pcross = m(2); % fraction cross peak (negative for diagonal peaks)
rdnp = m(3); % rate dnp buildup
rsd = m(4); % rate spin diffusion buildup
rcp = m(5); % rate cp buildup
r1 = m(6); % rate spin diffusion decay
r1rho = m(7); % rate cp decay

s = s0.*(1 - exp(-rdnp.*xps.d1))... %dnp
    .*(exp(-r1.*xps.d8) - pcross*exp(-rsd.*xps.d8))./(1-r1./rsd)... %spin diffusion
    .*(exp(-r1rho.*xps.p15) - exp(-rcp.*xps.p15))./(1-r1rho./rcp); %cp

