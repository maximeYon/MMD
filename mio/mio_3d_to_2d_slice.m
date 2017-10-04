function J = mio_3d_to_2d_slice(I, d, k, is_color)
% function J = mio_3d_to_2d_slice(I, d, k, is_color)

if (nargin < 2), d = 3; end
if (nargin < 3), k = []; end
if (nargin < 4), is_color = 0; end

if (is_color)
    I = permute(I, [2 3 4 1]);
end

if (isempty(k)), k = round(size(I, d)/2); end

switch (d)
    case 1
        J = I(k,:,:,:);
    case 2
        J = I(:,k,:,:);
    case 3
        J = I(:,:,k,:);
    otherwise
        error('check this');
end

J = squeeze(J);
J = J(:,end:-1:1,:);
J = permute(J, [2 1 3]);

