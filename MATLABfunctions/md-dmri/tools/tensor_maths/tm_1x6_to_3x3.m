function t = tm_1x6_to_3x3(t)
% function t = tm_1x6_to_3x3(t)
%
% convert a second order tensor from Voigt format to 3x3 matrix format


t = t .* [1 1 1 sqrt(1/2) sqrt(1/2) sqrt(1/2)];

t = t([1 4 5; 4 2 6; 5 6 3]);