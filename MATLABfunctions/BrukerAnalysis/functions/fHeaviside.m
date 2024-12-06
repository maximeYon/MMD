function Y = fHeaviside(X)
 
Y = zeros(size(X));
Y(X>0) = 1;