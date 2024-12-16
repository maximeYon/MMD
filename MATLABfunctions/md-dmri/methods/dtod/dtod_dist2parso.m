function [n,dparo,dperpo,theta,phi,w] = dtod_dist2parso(dtod,omega)

[n,dpar,dperp,theta,phi,d0,rpar,rperp,w] = dtod_dist2par(dtod);

asize = [1 numel(omega)];
dtodpars.omega = repmat(omega,[n 1]);
dtodpars.dpar = repmat(dpar,asize);
dtodpars.dperp = repmat(dperp,asize);
dtodpars.theta = repmat(theta,asize);
dtodpars.phi = repmat(phi,asize);
dtodpars.d0 = repmat(d0,asize);
dtodpars.rpar = repmat(rpar,asize);
dtodpars.rperp = repmat(rperp,asize);
dtodpars.w = repmat(w,asize);

dtodpars.dparo = dtodpars.d0 - (dtodpars.d0 - dtodpars.dpar)./(1 + dtodpars.omega.^2./dtodpars.rpar.^2);
dtodpars.dperpo = dtodpars.d0 - (dtodpars.d0 - dtodpars.dperp)./(1 + dtodpars.omega.^2./dtodpars.rperp.^2);

dparo = dtodpars.dparo;
dperpo = dtodpars.dperpo;
theta = dtodpars.theta;
phi = dtodpars.phi;
w = dtodpars.w;
