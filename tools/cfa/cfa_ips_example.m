function ips = cfa_ips_example()
% function ips = cfa_ips_example()

ips.B0 = 3;                             % main magnetic field [T]
ips.gm = msf_const_gamma;               % gyromagnetic constant [rad/T]

v_f = [1 0 0];                          % vector along freq. enc direction
v_p = [0 1 0];                          % vector along phase enc direction
v_s = [0 0 1];                          % vector along slice normal direction

ips.o   = [v_f; v_p; v_s];              % orientation matrix

ips.fov = [.2 .3 .1];                   % fov size along FPS [m]
ips.res = [2 2 4] * 1e-3;               % spatial resolution in FPS [m] 
ips.ipa = 2;                            % in-plane acceleration factor [1]
ips.ecs = .65e-3;                       % echo spacing [s]

ips.kpv = ips.ipa/ips.fov(2)/ips.ecs;   % k-phase velocity [1/(s.m)]

[X, Y, Z] = fov2xyz(ips.fov, ips.res); 
ips.r_xyz = [X(:) Y(:) Z(:)];           % position of all voxel centers [m]


% Not really imaging parameter -- but required for CFA
ips.T2s = 40e-3;                        % T2-star time of tissue [s]


