function msf_imagesc(I,d,k,c)
% function msf_imagesc(I,d,k,c)
%
% Displays a 2D slice through a 3D or 4D image volume I
%
% d - dimension
% k - slice
% c - volume

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


