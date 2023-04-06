function I = mio_min_max_cut(I, min_max, tmp)
% function I = mio_min_max_cut(I, min_max)
%
% cut values below min or above max, where [min max] = min_max

if (nargin == 3), min_max = [min_max tmp]; end

I = min(max(I, min_max(1)), min_max(2));
