function I = mio_min_max_cut(I, min_max)
% function I = mio_min_max_cut(I, min_max)
%
% cut values below min or above max, where [min max] = min_max

I = min(max(I, min_max(1)), min_max(2));
