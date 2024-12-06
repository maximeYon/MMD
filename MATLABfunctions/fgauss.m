function P = fgauss(x,meanx,stdx)

P = 1/stdx/sqrt(2*pi)*exp(-(x-meanx).^2/2/stdx^2);