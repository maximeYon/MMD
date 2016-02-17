function M = mio_mask_keep_largest(M)
% function M = mio_mask_keep_largest(M)

M = M > 0;

if (all(M(:) == 0))
    warning('mask is empty');
    return; 
end

% Remove isolated areas
[L,n] = bwlabeln(M > 0,6);

n_elem = zeros(1,n);
for c = 1:n
    n_elem(c) = sum(L(:) == c);
end

[~,ind_max] = max(n_elem);

M = L == ind_max;





