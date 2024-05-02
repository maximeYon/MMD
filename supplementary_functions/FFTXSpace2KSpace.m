function [PostFFT] = FFTXSpace2KSpace(PreFFT, Dim)
   
   PostFFT = fftshift(ifft(ifftshift(PreFFT, Dim), [], Dim), Dim) ;
%    PostFFT = ifft(PreFFT, [], Dim) ;

return ;