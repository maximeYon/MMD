function [dparo,dperpo,theta,phi,r1,r2,w] = dtor1r2d_4d_m2parso(m,omega)

[dpar,dperp,theta,phi,d0,rpar,rperp,r1,r2,w] = dtor1r2d_4d_m2pars(m);

dparo = d0 - (d0 - dpar)./(1 + omega.^2./rpar.^2);
dperpo = d0 - (d0 - dperp)./(1 + omega.^2./rperp.^2);
