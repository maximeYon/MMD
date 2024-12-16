function iso = tm_eigvals2iso(eigvals_nx3)
% function iso = tm_eigvals2iso(eigvals_nx3)
%
% Input: nx3 matrix of eigenvalues
% Output: nx1 vector of isotropic average values

eigval1 = eigvals_nx3(:,1);
eigval2 = eigvals_nx3(:,2);
eigval3 = eigvals_nx3(:,3);

iso = (eigval1+eigval2+eigval3)/3;
