function [M, D, dD, A] = ilt_regularized_inversion_kernel(xps)

% Extract b-values
b = xps.b;

% Sample diffusivities to build the inversion kernel matrix
M = 150;
D = linspace(0,3e-9,M+1)';
dD = diff(D);
D = D(1:M) + dD/2;

% Build inversion kernel matrix
[~,Darray] = ndgrid(b,D);
[barray,dDarray] = ndgrid(b,dD);
A = 1./Darray.*dDarray.*exp(-barray.*Darray);
