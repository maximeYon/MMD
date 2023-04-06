function fa = tm_eigvals2fa(eigvals_nx3)
% function fa = tm_eigvals2fa(eigvals_nx3)
%
% Input: nx3 matrix of eigenvalues
% Output: nx1 vector of FA values

eigval1 = eigvals_nx3(:,1);
eigval2 = eigvals_nx3(:,2);
eigval3 = eigvals_nx3(:,3);

iso = tm_eigvals2iso(eigvals_nx3);
fa = sqrt((3/2)*((eigval1-iso).^2+(eigval2-iso).^2+(eigval3-iso).^2)...
    ./(eigval1.^2+eigval2.^2+eigval3.^2));