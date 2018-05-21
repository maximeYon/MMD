function ips = cfa_ips_update_derived_pars(ips)

[X, Y, Z] = fov2xyz(ips.fov, ips.res); 
ips.r_xyz = [X(:) Y(:) Z(:)];           % Location of voxels

ips.kpv = ips.ipa/ips.fov(2)/ips.ecs;   % k-phase velocity [1/(s.m)]