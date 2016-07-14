function msf_imagesc(I,d,k,c)
% function msf_imagesc(I,d,k,c)
%
% Displays a 2D slice through a 3D or 4D image volume I
%
% d - dimension
% k - slice
% c - volume

if (size(I,1) == 3) && (ndims(I) == 4)
    % assume it is a color image
    
    if (nargin < 2), d = 3; end
    if (nargin < 3) || isempty(k), k = round(size(I, d + 1)/2 ); end
    
    switch (d)
        case 1
            tmp = I(:,k,:,:);
        case 2
            tmp = I(:,:,k,:);
        case 3
            tmp = I(:,:,:,k);
    end
    
    
    if (max(tmp(:)) > 1)
        tmp = tmp / quantile(tmp(:), 0.99);
    end
    tmp(tmp > 1) = 1;
    tmp(tmp < 0) = 0;
    
    imagesc(permute(squeeze(tmp(:,:,end:-1:1)), [3 2 1]));
    axis image off;
            
    
else
    
    if (nargin < 2), d = 3; end
    if (nargin < 3) || (isempty(k)), k = round(size(I,d)/2); end
    if (nargin < 4) || (isempty(c)), c = 1; end
    
    colormap(gray);
    
    switch (d)
        case 1
            tmp = I(k,:,:,c);
        case 2
            tmp = I(:,k,:,c);
        case 3
            tmp = I(:,:,k,c);
    end
    imagesc(flipud(squeeze(mean(tmp,d))'));
    
    axis image off;
    
    
end