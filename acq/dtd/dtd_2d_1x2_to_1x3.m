function t = dtd_2d_1x2_to_1x3(ad, rd, n)
% function t = dtd_2d_1x2_to_1x3(ad, rd, n)
%
% Two-dimensional version of dtd_1x3_to_1x6

x = n(:,1);
y = n(:,2);

o = zeros(size(x));
e = ones(size(x));
c = sqrt(2);

t = [x.*x y.*y c*x.*y] * (ad - rd) + ...
    [e e o ] * rd;
