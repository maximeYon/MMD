function Y = fMurdayCotts(Pin,Xin,Pnorm,Xnorm,Ynorm,alphaR,D,delta,DELTA);Pin = Pin.*Pnorm;G2 = Xin'*Xnorm;Y0 = Pin(1);R = Pin(2);sigma = Pin(3);Rvector = linspace(R-3*sigma,R+3*sigma);dR = Rvector(2) - Rvector(1);gamma = 26.75e7;[alphaRarray,G2array,Rarray] = ndgrid(alphaR,G2,Rvector);alphaarray = alphaRarray./Rarray;PR = 1/sigma/sqrt(2*pi)*exp(-0.5*((Rarray-R)/sigma).^2);YR = Y0*exp(sum(-2*gamma^2*G2array/D.*(alphaarray.^-4./(alphaarray.^2.*Rarray.^2-2).*(2*delta-(2+exp(-alphaarray.^2*D*(DELTA-delta))-2*exp(-alphaarray.^2*D*delta)-2*exp(-alphaarray.^2*D*DELTA)+exp(-alphaarray.^2*D*(DELTA+delta)))./alphaarray.^2/D))));Y = sum(PR(1,:,:).*YR*dR,3);Y = Y'/Ynorm;