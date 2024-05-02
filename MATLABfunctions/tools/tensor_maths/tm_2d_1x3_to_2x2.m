function t = tm_2d_1x3_to_2x2(t)
% function t = tm_2d_1x3_to_2x2(t)
%
% Two-dimensional version of tm_1x6_to_3x3


t = t([1 3;3 2]) ./ sqrt([1 2;2 1]);