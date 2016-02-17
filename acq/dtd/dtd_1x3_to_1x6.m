function t = ts_make_tensors_from_directions(ad, rd, n)
% function DT_dist = ts_make_tensors_from_directions(ad, rd, n)

x = n(:,1);
y = n(:,2);
z = n(:,3);

o = zeros(size(x));
e = ones(size(x));
c = sqrt(2);

if (numel(ad) * numel(rd) == 1)
    
    t = [x.*x y.*y z.*z c*x.*y c*x.*z c*y.*z] * (ad - rd) + ...
        [e e e o o o ] * rd;
    
else
    
    t = [x.*x y.*y z.*z c*x.*y c*x.*z c*y.*z] .* repmat(ad - rd,[1 6]) + ...
        [e e e o o o ] .* repmat(rd, [1 6]);
    
    
end
