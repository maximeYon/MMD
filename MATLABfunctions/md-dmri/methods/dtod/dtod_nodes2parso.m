function [n,dparo,dperpo,theta,phi] = dtod_nodes2parso(dtod_nodes,omega)

[n,dpar,dperp,theta,phi,d0,rpar,rperp] = dtod_nodes2par(dtod_nodes);

asize = [1 numel(omega)];
dtodpars.omega = repmat(omega,[n 1]);
dtodpars.dpar = repmat(dpar,asize);
dtodpars.dperp = repmat(dperp,asize);
dtodpars.theta = repmat(theta,asize);
dtodpars.phi = repmat(phi,asize);
dtodpars.d0 = repmat(d0,asize);
dtodpars.rpar = repmat(rpar,asize);
dtodpars.rperp = repmat(rperp,asize);

dtodpars.dparo = dtodpars.d0 - (dtodpars.d0 - dtodpars.dpar)./(1 + dtodpars.omega.^2./dtodpars.rpar.^2);
dtodpars.dperpo = dtodpars.d0 - (dtodpars.d0 - dtodpars.dperp)./(1 + dtodpars.omega.^2./dtodpars.rperp.^2);

dparo = dtodpars.dparo;
dperpo = dtodpars.dperpo;
theta = dtodpars.theta;
phi = dtodpars.phi;
