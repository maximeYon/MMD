function [ir,jr,kr] = mio_mask_ranges(M,opt)
% function [ir,jr,kr] = mio_mask_ranges(M,opt)
%
% Detect the "bounding box" of a mask

ind = @(d1,d2) squeeze(sum( sum( M, d1), d2) > 0);

ir = find( ind(2,3), 1, 'first'):find( ind(2,3), 1, 'last');
jr = find( ind(1,3), 1, 'first'):find( ind(1,3), 1, 'last');
kr = find( ind(1,2), 1, 'first'):find( ind(1,2), 1, 'last');


