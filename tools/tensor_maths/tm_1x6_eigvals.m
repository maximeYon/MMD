function lambda = tm_1x6_eigvals(G)
% function lambda = tm_1x6_eigvals(G)
%
% xx yy zz yz xz xy
% 
% xx yy zz xy xz yz

% Use Cardano's method to diagonalization
c_0 = ...
    G(:,1) .* abs(G(:,6)).^2/2 + ...
    G(:,2) .* abs(G(:,5)).^2/2 + ...
    G(:,3) .* abs(G(:,4)).^2/2 - ...
    G(:,1) .* G(:,2) .* G(:,3) - ...
    2 * real(conj(G(:,5)) .* G(:,6) .* G(:,4)) / sqrt(2)^3;

c_1 = G(:,1).*G(:,2) + G(:,1).*G(:,3) + G(:,2).*G(:,3) - ...
    sum(abs(G(:,4:6)).^2/2,2);

c_2 = -sum(G(:,1:3),2);


p = c_2.^2 - 3 * c_1;
q = -27/2 * c_0 - c_2.^3  + 9/2 * c_2 .* c_1;

phi = 1/3 * atan( sqrt( 27 * (1/4 * c_1.^2 .* (p - c_1) + c_0 .* (q + 27/4 * c_0))) ./ q);

phi = phi + pi * (q < 0); 


x1 = 2 * cos(phi);
x2 = -cos(phi) - sqrt(3) * sin(phi);
x3 = -cos(phi) + sqrt(3) * sin(phi);

ps = sqrt(p) / 3;

lambda = [ps .* x1 - 1/3 * c_2 ps .* x2 - 1/3 * c_2 ps .* x3 - 1/3 * c_2] ;




