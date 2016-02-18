function t = dtd_1x6_to_3x3(t)
% function t = dtd_1x6_to_3x3(t)
%
% convert a second order tensor from Voigt format to 3x3 matrix format

c = sqrt(2);
t = t .* [1 1 1 c c c].^(-1);

t = t([1 4 5;4 2 6;5 6 3]);