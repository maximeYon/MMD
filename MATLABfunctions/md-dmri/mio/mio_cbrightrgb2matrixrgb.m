function I = mio_cbrightrgb2matrixrgb(c)
% function I = mio_cbrightrgb2matrixrgb(c)
%
% Converts structure c with fields 'bright', 'r', 'g', and 'b'
% to image matrix.
% Replaces NaN and Inf with 0 and cuts intensities outside the range 0-1

I = zeros([size(c.bright) 3]);
I(:,:,1) = c.bright .* c.r;
I(:,:,2) = c.bright .* c.g;
I(:,:,3) = c.bright .* c.b;

I(isnan(I)) = 0;
I(isinf(I)) = 0;

I( I(:) > 1 ) = 1;
I( I(:) < 0 ) = 0;

