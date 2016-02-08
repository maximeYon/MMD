function t = ts_2d_1x3_to_2x2(t)

t = t([1 3;3 2]) ./ sqrt([1 2;2 1]);