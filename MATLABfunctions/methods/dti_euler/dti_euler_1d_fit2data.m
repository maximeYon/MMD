function s = dti_euler_1d_fit2data(m, xps)
% function s = dti_euler_1d_fit2data(m, xps)


% convert to readable parameters
s0          = m(1);
lambda_x    = m(2);
lambda_y    = m(3);
lambda_z    = m(4);
euler_alpha = m(5);
euler_beta  = m(6);
euler_gamma = m(7);

% from euler angles to rotation matrices
[rotmat,rotmatinv] = tm_euler_angles2rotmat(euler_alpha,euler_beta,euler_gamma);

% dt in principal axis system (PAS)
dt_pas = [
    lambda_x 0 0
    0 lambda_y 0
    0 0 lambda_z];

% active rotation of dt from PAS = lab
dt = rotmat*dt_pas*rotmatinv;

% convert dt to voigt
dt = tm_3x3_to_1x6(dt);

% signal for all b-tensors
s = m(1) * exp(-(dt * xps.bt'))';
