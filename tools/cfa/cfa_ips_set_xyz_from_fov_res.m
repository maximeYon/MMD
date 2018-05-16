function ips = cfa_ips_set_xyz_from_fov_res(ips)
% function ips = cfa_ips_set_xyz_from_fov_res(ips)

[X, Y, Z] = fov2xyz(ips.fov, ips.res); 
ips.r_xyz = [X(:) Y(:) Z(:)];