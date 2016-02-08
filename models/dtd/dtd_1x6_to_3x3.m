function t = ts_1x6_to_3x3(t)

c = sqrt(2);
t = t .* [1 1 1 c c c].^(-1);

t = t([1 4 5;4 2 6;5 6 3]);