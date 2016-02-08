function t = ts_3x3_to_1x6(t)

c = sqrt(2);
t = t([1 5 9 2 3 6]) .* [1 1 1 c c c];