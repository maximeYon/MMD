function I = fLorentz(omega,tau)

I = 2*tau./(1+omega.^2.*tau.^2);