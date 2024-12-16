function [dpar,dperp,theta,phi,d0,rpar,rperp,w] = dtod_bsmerge_catpars(dpar,dperp,theta,phi,d0,rpar,rperp,w,dpar_temp,dperp_temp,theta_temp,phi_temp,d0_temp,rpar_temp,rperp_temp,w_temp)

dpar = cat(4,dpar,dpar_temp);
dperp = cat(4,dperp,dperp_temp);
theta = cat(4,theta,theta_temp);
phi = cat(4,phi,phi_temp);
d0 = cat(4,d0,d0_temp);
rpar = cat(4,rpar,rpar_temp);
rperp = cat(4,rperp,rperp_temp);
w = cat(4,w,w_temp);
