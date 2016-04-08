function [rotmat,rotmatinv] = tm_euler_angles2rotmat(alpha,beta,gamma)
% function [rotmat,rotmatinv] = tm_euler_angles2rotmat(alpha,beta,gamma)
%
% Active right-handed rotation
% First gamma around z, then beta around y, and finally alpha around z
% Input must be scalars.


rotmat_alpha = [
    cos(alpha) -sin(alpha) 0
    sin(alpha) cos(alpha) 0
    0           0          1];

rotmat_beta = [
    cos(beta)  0 sin(beta)
    0          1 0
    -sin(beta) 0 cos(beta)];

rotmat_gamma = [
    cos(gamma)  -sin(gamma) 0
    sin(gamma)  cos(gamma) 0
    0           0          1];

rotmat = rotmat_alpha*rotmat_beta*rotmat_gamma;
rotmatinv = inv(rotmat);