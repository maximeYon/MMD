
function epsilon1 = findzeroPetrache(msx, epsilon1guess)

options = optimset('Display', 'off','TolX',eps);
epsilon1 = fzero(@Eq30, epsilon1guess, options);

 function y = Eq30(epsilon1)
 y = 1 + 2./epsilon1.^2 - 2./epsilon1.*coth(epsilon1) - msx;
 end
end
