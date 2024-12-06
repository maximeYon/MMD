
function y = findBesselzeroGPD(nD, x0)

options = optimset('Display', 'off','TolX',eps);
y = fzero(@bessels, x0, options);

 function y = bessels(x) % Compute the polynomial.
 y = besselj(nD/2,x) - x.*besselj(nD/2+1,x);
 end
end
