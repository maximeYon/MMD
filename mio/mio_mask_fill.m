function M = mio_mask_fill(M,d)
% function M = mio_mask_fill(M,d)
%
% fill holes in a mask on a slice-by-slice basis

if (nargin < 2), d = 3; end

for k = 1:size(M,d)
    
    switch (d)
        case 1
            M_tmp = squeeze(M(k,:,:));
        case 2
            M_tmp = squeeze(M(:,k,:));
        case 3
            M_tmp = squeeze(M(:,:,k));
    end
    
    % Create a larger map
    I_tmp = ones(size(M_tmp) + 2);
    I_tmp(2:(end-1),2:(end-1)) = 1 - M_tmp;
    [L,n] = bwlabel(I_tmp);
    L = L(2:(end-1), 2:(end-1));
    
    % Continue
    n_elem = zeros(1,n);
    for c = 1:n, n_elem(c) = sum(L(:) == c); end
    for c = 1:n
        if (n_elem(c) ~= max(n_elem))
            M_tmp(L == c) = 1;
        end
    end
    
    switch (d)
        case 1
            M(k,:,:) = M_tmp;
        case 2
            M(:,k,:) = M_tmp;
        case 3
            M(:,:,k) = M_tmp;
    end
end



