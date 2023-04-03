function ans=erfi(x)
% %erfi(x). The Imaginary error function, as it is defined in Mathematica
% %erfi(z)==erf(iz)/i (z could be complex) using 
% %the incomplete gamma function in matlab: gammainc
% %Using "@": erfi = @(x) real(-sqrt(-1).*sign(x).*gammainc(-x.^2,1/2))
% %Note: limit(x->0) erfi(x)/x -> 2/sqrt(pi)
%
% %Example 1: 
% x=linspace(0.001,6,100);
% y=exp(-x.^2).*erfi(x)./2./x;
% figure(1), clf;plot(x,y*sqrt(pi))
%
% %Example 2: 
% [x,y]=meshgrid(linspace(-3,3,180),linspace(-3,3,180));
% z=x+i*y;
% figure(1), clf;contourf(x,y,log(erfi(z)))
% axis equal;axis off
xc=5.7;%cut for asymptotic approximation (when x is real)
ans=~isreal(x).*(-(sqrt(-x.^2)./(x+isreal(x))).*gammainc(-x.^2,1/2))+...
    isreal(x).*real(-sqrt(-1).*sign(x).*((x<xc).*gammainc(-x.^2,1/2))+...
    (x>=xc).*exp(x.^2)./x/sqrt(pi));

% function y=erfi(z)
% % function y=erfi( imag(z) ) = imag(erf(z))
% % Returns imaginary part of the error function of imaginary z
% 
% % ERFI.M 4/30/92	 Hilmar Schlegel 
% % Internet: hshlgaii@mailszrz.zrz.tu-berlin.de (52.32/13.25)
% % Latest change: 5/6/92
% 
% t=real(z); % arg check
% if imag(z),disp('*** Warning: Imaginary part not real in ERFI '),end 
% [nz,mz]=size(z);
% t=t(:);
% s1=t;
% s2=t;
% z=t.^2;
% new=(1:length(t))';
% n=1;
% while 1
% t(new)=t(new).*z(new)/n;
% s2(new)=s1(new)+t(new)/(2*n+1);
% new=find(s2~=s1);
% if isempty(new),break,end
% n=n+1;
% s1(new)=s2(new);
% end
% y=reshape(2*s2/sqrt(pi),nz,mz);