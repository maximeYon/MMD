% Gradient waveforms for axisymmetric b-tensors with 
% variable frequency content by double rotation of the q-vector.
%
% Work in progress
%

% Define timing parameters relative to a total echo time = 1
np = 1000; % Number of time intervals in waveform [1000]
% Double cone parameters
% q-vector buildup parameters

out_path = '/opt/topspin4.0.7/exp/stan/nmr/lists/gp/user';
dorratios = 0:5;

epsilon_up = .03; plateau = 0; epsilon_down = 4*epsilon_up;

% Parameters for calculating b-values
gmax = 100/100*3; % Max gradient [3 T/m]
tau = 25e-3; % Waveform duration [10 ms]

%------------------------

% Define path for output folders
%[run_path,run_name,run_ext] = fileparts(run_fn);
run_path = pwd;
%out_path = fullfile(run_path,['protocol_dor']);

for nshape = 1:length(dorratios)
    dorratio = dorratios(nshape);
    

    if dorratio == 0
        zeta1 = 54.74/180*pi; zeta2 = 0*90/180*pi;
    else
        zeta1 = -54.74/180*pi; zeta2 = 1*90/180*pi;
    end
    deltapsi1 = 1*2*pi; deltapsi2 = dorratio*2*pi;

    taured = 1;
    dt = taured/np;
    t = taured*linspace(0,1,np);

    % Quarter-sine ramp up
    np_epsilon_up = round(epsilon_up/dt);
    t_epsilon_up = pi/2*linspace(0,1,np_epsilon_up)';
    g_up = sin(t_epsilon_up);
    %figure(1), clf, plot(t_epsilon_up,g_up,'-'), return

    % Half-sine ramp down
    np_epsilon_down = round(epsilon_down/dt);
    t_epsilon_down = pi/2*linspace(-1,1,np_epsilon_down)';
    g_down = 1 - .5*(1+sin(t_epsilon_down));
    %figure(1), clf, plot(t_epsilon_down,g_down,'-'), return

    % Plateau
    np_plateau = round(plateau/dt);
    np_inter = round(taured/dt)-2*(np_epsilon_up+np_plateau+np_epsilon_down);

    ga = [g_up; ones(np_plateau,1); g_down; zeros(np_inter,1); g_down-1; -1*ones(np_plateau,1); -flipud(g_up)];
    %figure(1), clf, plot(t,ga,'-'), return

    % Some equations taken from
    % Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016).
    % http://dx.doi.org/10.1039/c5cp07251d

    q = cumsum(ga); % Topgaard Eq 35 (Note: all dt cancel and can be omitted.)
    q = mean([q flipud(q)],2); % Removes asymmetry introduced by cumsum
    q = q/max(q);
    b = sum(q.^2)*dt; % Topgaard Eq 32
    psi_norm = 1/b*cumsum(q.^2)*dt; % Topgaard Eq 31
    psi_norm = mean([psi_norm 1-flipud(psi_norm)],2); % Removes asymmetry introduced by cumsum

    psi0 = 0;
    psi = psi0 + deltapsi1*psi_norm; % Topgaard Eq 31

    qx = zeros(np,1);
    qy = zeros(np,1);
    qz = zeros(np,1);

    q0 = [0; 0; 1];

    % Double rotation of q-vector by two sets of Euler angles alpha, beta, gamma
    for count = 1:np
        alpha2 = deltapsi2/deltapsi1*psi(count);
        beta2 = zeta2;
        gamma2 = 0;
        [rotmat2,rotmatinv2] = tm_euler_angles2rotmat(alpha2,beta2,gamma2);
        alpha1 = psi(count);
        beta1 = zeta1;
        gamma1 = 0;
        [rotmat1,rotmatinv1] = tm_euler_angles2rotmat(alpha1,beta1,gamma1);
        q_temp = q(count)*rotmat1*rotmat2*q0;
        %q_temp = rotmat2*rotmat1*q0;
        qx(count) = q_temp(1);
        qy(count) = q_temp(2);
        qz(count) = q_temp(3);
    end
    % figure(1), clf, plot(1:np,qx,'r-'), return

    % Assure that first and last values equal zero
    qx([1 end]) = 0;
    qy([1 end]) = 0;
    qz([1 end]) = 0;

    gx = gradient(qx)/dt;
    gy = gradient(qy)/dt;
    gz = gradient(qz)/dt;

    % Assure that first and last values equal zero
    gx([1 end]) = 0;
    gy([1 end]) = 0;
    gz([1 end]) = 0;

    re_gr = gx*sqrt(3/2);
    im_gr = gy*sqrt(3/2);
    ga = gz*sqrt(3);

    gnorm = max(abs([ga(:); re_gr(:)+1i*im_gr(:); ...
        sqrt(re_gr(:).^2*2/3+im_gr(:).^2*2/3+ga(:).^2*1/3)]));

    ga = ga/gnorm;
    re_gr = re_gr/gnorm;
    im_gr = im_gr/gnorm;
    gx = gx/gnorm;
    gy = gy/gnorm;
    gz = gz/gnorm;
    gr = sqrt(gx.^2+gy.^2+gz.^2);

    qx = cumsum(gx);
    qy = cumsum(gy);
    qz = cumsum(gz);
    bxx = sum(qx.^2);
    byy = sum(qy.^2);
    bzz = sum(qz.^2);
    bxy = sum(qx.*qy);
    bxz = sum(qx.*qz);
    byz = sum(qy.*qz);

    b_tensor = [[bxx bxy bxz]; [bxy byy byz]; [bxz byz bzz]]/(bxx+byy+bzz)
    % Should be diagonal and all eigenvalues 1/3

    gamma = 26.75e7;
    dt = tau/np;
    t = tau*linspace(0,1,np)';
    q = gamma*gmax*cumsum(ga*dt);
    b = sum(q.^2*dt);
    qmax = max(q);
    td = b/qmax^2;
    C = b/gamma^2/gmax^2/tau^3;
    Cref = 3.1966e-4; % for dorratio = 5
    cnst1 = sqrt(Cref/C)*100;

    dgadt = gradient(gmax*ga/dt);
    dgrdt = gradient(gmax*abs(gr)/dt);
    dre_grdt = gradient(gmax*re_gr/dt);
    dim_grdt = gradient(gmax*im_gr/dt);

    figure(1), clf
    subplot(2,1,1)
    plot(t,gmax*re_gr,'r-',t,gmax*im_gr,'g-',t,gmax*ga,'b-',t,gmax*abs(gr),'k--')
    ylabel('g / Tm^-^1')
    title(['b = ' num2str(b/1e9,2) '\cdot10^9 sm^-^2   q = ' num2str(qmax/2/pi/1e6,2) '\cdot10^6 m^-^1'])

    subplot(2,1,2)
    plot(t,dre_grdt,'r-',t,dim_grdt,'g-',t,dgadt,'b-',t,dgrdt,'k--')
    xlabel('t / s'), ylabel('(dg/dt) / Tm^-^1s^-^1')

    res = bruker_mkshapefile(fullfile(out_path,['ga' num2str(nshape)]),re_gr*cnst1/100,['dor' num2str(dorratio) ' radial a']);
    res = bruker_mkshapefile(fullfile(out_path,['gb' num2str(nshape)]),im_gr*cnst1/100,['dor' num2str(dorratio) ' radial b']);
    res = bruker_mkshapefile(fullfile(out_path,['gc' num2str(nshape)]),ga*cnst1/100,['dor' num2str(dorratio) ' axial c']);
    
    pause(1)
end
