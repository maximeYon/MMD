clear

% Here are some examples of asymmetric diffusion encoding in a spin-echo
% sequence. See the example cases to experiment with some common issues.

example_nr = 5;

switch example_nr
    case 1
        % The stejskal-tanner experiment is symmetric and will never
        % produce any maxwell terms, regardless of rotations of gwf and
        % fov.
        wf_type    = 1; % ST
        find_worst = 2; % none of the modes have any effect
        
    case 2
        % If we run a DDE with orthogonal gradients in the xy-plane, no
        % errors will appear!
        wf_type    = 2; % DDE in x and y
        find_worst = 3; % rotating the fov has no effect
        
    case 3
        % If we use the DDE in x and y, but allow it to rotate, there will
        % be a gross signal error away from the iso center.
        wf_type    = 2; % DDE in x and y
        find_worst = 2; % worst rotation of gradients
        
    case 4
        % If we start with DDE in the xz-plane, there will always be an
        % error, even if the gradient and the fov are not rotated.
        wf_type    = 3; % DDE in x and z
        find_worst = 1; % no rotation
        
    case 5
        % If we start with DDE in the xz-plane, there will always be an
        % error, and it can get even worse if we search for the worst fov
        % or gvf rotation.
        wf_type    = 3; % DDE in x and z
        find_worst = 3; % no rotation
end


% Fetch imaging parameter structure (ips), and gradient waveform (gwf)
ips = cfa_ips_example();
[gwf, rf, dt] = cfa_gwf_example(wf_type); % 1 is ST-SDE; 2 is ortho DDE in x & y; 3 is ortho DDE in x & z

% Rotations of the FOV and the gwf are imortant and we can fing the worst
% rotation by an fminsearch.
switch find_worst
    case 1 % let it be
        
    case 2 % Find worst rotation of gwf
        R = cfa_find_worst_gwf_rot(gwf, rf, dt, ips, 1);
        gwf = gwf * R;
        
    case 3 % Find worst rotation of FOV
        R = cfa_find_worst_fov_rot(gwf, rf, dt, ips, 1);
        ips = cfa_apply_R_to_fov(ips, R);
        
end

% Calculate the bias field in the volume defined by ips.
c = cfa_maxwell_bias(gwf, rf, dt, ips);

% We can also calculate the bias field for arbitrary points in space by
% switching out ips.r_xyz for custom coordinates. For example, we can
% calculate the bias on a convex surface (here an ellipsoid).
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

