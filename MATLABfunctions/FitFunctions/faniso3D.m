function Y = faniso3D(Pin,Xin,Pnorm,Xnorm,Ynorm);if nargin == 2    Pnorm = ones(size(Pin));    Xnorm = 1;    Ynorm = 1;endPin = Pin.*Pnorm;Xin = Xin*Xnorm;I0 = Pin(1);Dpar = Pin(2);Dperp = Pin(3);b = Xin;if Dpar>Dperp    E = sqrt(pi)/2*exp(-b*Dperp)./sqrt(b*(Dpar-Dperp)).*erf(sqrt(b*(Dpar-Dperp)));elseif Dperp>Dpar    E = sqrt(pi)/2*exp(-b*Dperp)./sqrt(b*(Dpar-Dperp))*i.*erfi(imag(sqrt(b*(Dpar-Dperp))));else    E = exp(-b*Dperp);endY = I0*E;% N = 100;% x = linspace(0,1,N);% % [karray, xarray] = ndgrid(Xin,x);% % SD = exp(-karray.*Dperp).*exp(-karray.*(Dpar-Dperp).*xarray.^2);% % mSD = (SD(:,1:(N-1)) + SD(:,2:N))/2;% dx = x(2:N) - x(1:(N-1));% % Y = I0.*mSD*dx';Y = Y/Ynorm;