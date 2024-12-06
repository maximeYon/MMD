function [X,Y] = fSchmidt(x,y,z)

%Check normalization
r = sqrt(x.^2 + y.^2 + z.^2);

if any(any(abs(r-1)>.001))
    warning('Input vector not normalized')
end

X = -sqrt(2./(1+z)).*x;
Y = -sqrt(2./(1+z)).*y;
