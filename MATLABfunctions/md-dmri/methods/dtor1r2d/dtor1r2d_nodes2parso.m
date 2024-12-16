function [n,dparo,dperpo,theta,phi,r1,r2] = dtor1r2d_nodes2parso(dtor1r2d_nodes,omega)

[n,dpar,dperp,theta,phi,d0,rpar,rperp,r1,r2] = dtor1r2d_nodes2par(dtor1r2d_nodes);

asize = [1 numel(omega)];
dtor1r2dpars.omega = repmat(omega,[n 1]);
dtor1r2dpars.dpar = repmat(dpar,asize);
dtor1r2dpars.dperp = repmat(dperp,asize);
dtor1r2dpars.theta = repmat(theta,asize);
dtor1r2dpars.phi = repmat(phi,asize);
dtor1r2dpars.d0 = repmat(d0,asize);
dtor1r2dpars.rpar = repmat(rpar,asize);
dtor1r2dpars.rperp = repmat(rperp,asize);
dtor1r2dpars.r1 = repmat(r1,asize);
dtor1r2dpars.r2 = repmat(r2,asize);

dtor1r2dpars.dparo = dtor1r2dpars.d0 - (dtor1r2dpars.d0 - dtor1r2dpars.dpar)./(1 + dtor1r2dpars.omega.^2./dtor1r2dpars.rpar.^2);
dtor1r2dpars.dperpo = dtor1r2dpars.d0 - (dtor1r2dpars.d0 - dtor1r2dpars.dperp)./(1 + dtor1r2dpars.omega.^2./dtor1r2dpars.rperp.^2);

dparo = dtor1r2dpars.dparo;
dperpo = dtor1r2dpars.dperpo;
theta = dtor1r2dpars.theta;
phi = dtor1r2dpars.phi;
r1 = dtor1r2dpars.r1;
r2 = dtor1r2dpars.r2;
