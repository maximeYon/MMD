function msf_imagesc(I,d,k,c)
% function msf_imagesc(I,d,k,c)
%
% Displays a 2D slice through a 3D or 4D image volume I
%
% d - dimension
% k - slice
% c - volume

if (nargin < 2) || (isempty(d)), d = 3; end
if (nargin < 3) || (isempty(k)), k = round(size(I,d)/2); end
if (nargin < 4) || (isempty(c)), c = 1; end

if (size(I,1) == 3) && (ndims(I) == 4) % assume it is a color image
    
    tmp = mio_3d_to_2d_slice(I,d,k,1);
    
    if (isa(tmp, 'uint8')), tmp = double(tmp) / 255; end
    
    if (max(tmp(:)) > 1)
        tmp = tmp / quantile(tmp(:), 0.999);
    end
    
    tmp(tmp > 1) = 1;
    tmp(tmp < 0) = 0;
    
else
    colormap(gray);
    tmp = mio_3d_to_2d_slice(I(:,:,:,c),d,k);
end

imagesc(tmp);
axis image off;
