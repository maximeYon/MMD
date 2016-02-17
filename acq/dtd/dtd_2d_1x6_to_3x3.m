function t = ts_2d_1x6_to_3x3(t)

t = t([...
    1 3 4;
    3 2 5;
    4 5 6]) .* sqrt(1./[1 2 2; 2 1 2;2 2 1]);