function c = cfa_gwf_to_error_map(gwf, rf, dt, ips, find_mode)
% function c = cfa_gwf_to_error_map(gwf, rf, dt, ips, find_mode)

if nargin < 4
    ips = cfa_ips_example();
end

if nargin < 5
    find_mode = 0;
end

% Rotations of the FOV and the gwf are imortant and we can fing the worst
% rotation by an fminsearch.
switch find_mode
    case 0 % let it be
        
    case 1 % Find worst rotation of gwf
        R = cfa_find_worst_gwf_rot(gwf, rf, dt, ips, 1);
        gwf = gwf * R;
        
    case 2 % Find worst rotation of FOV
        R = cfa_find_worst_fov_rot(gwf, rf, dt, ips, 1);
        ips = cfa_apply_R_to_fov(ips, R);
        
end

% Calculate the bias field in the volume defined by ips.
c = cfa_maxwell_bias(gwf, rf, dt, ips);

% We can also calculate the bias field for arbitrary points in space by
% switching out ips.r_xyz for custom coordinates. For example, we can
% calculate the bias on a surface (here an ellipsoid).
ips2 = ips;
ips2.r_xyz = cfa_ellipsoid_xyz_from_fov(ips);
cs = cfa_maxwell_bias(gwf, rf, dt, ips2);


%% PLOT RELATIVE SIGNAL
clf
subplot(2,1,1)
h1 = cfa_plot_bias_volume (c, ips);

subplot(2,1,2)
h2 = cfa_plot_bias_surface(cs, ips2);

linkprop([h1, h2], {'CameraPosition','CameraUpVector'})
colormap (flip(hot))
set(gcf, 'color', 'w')

