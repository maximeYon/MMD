function t = tm_1x3_to_1x6(ad, rd, n)
% function t = tm_1x3_to_1x6(ad, rd, n)
%
% Convert a vector (1x3) to a tensor (1x6), but enable some shaping of the
% tensor, so that its longest eigenvalue is 'ad' and its shortest to 
% eigenvalues are 'rd'

x = n(:,1);
y = n(:,2);
z = n(:,3);

o = zeros(size(x));
e = ones(size(x));
c = sqrt(2);

if (numel(ad) * numel(rd) == 1)
    
    t = [x.*x y.*y z.*z c*x.*y c*x.*z c*y.*z] * (ad - rd) + ...
        [e e e o o o ] * rd;
    
elseif (size(n,1) == 1)

    t = repmat([x.*x y.*y z.*z c*x.*y c*x.*z c*y.*z], numel(ad), 1) .* repmat(ad - rd, 1, 6) + ...
        repmat([e e e o o o ], numel(rd), 1) .* repmat(rd, 1, 6);
        
else
    
    t = [x.*x y.*y z.*z c*x.*y c*x.*z c*y.*z] .* repmat(ad - rd,[1 6]) + ...
        [e e e o o o ] .* repmat(rd, [1 6]);
        
end
