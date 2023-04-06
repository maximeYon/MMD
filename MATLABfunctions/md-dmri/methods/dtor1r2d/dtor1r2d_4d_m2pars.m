function [dpar,dperp,theta,phi,d0,rpar,rperp,r1,r2,w] = dtor1r2d_4d_m2pars(m)

sz = size(m);

ind = false(sz(4),1);
ind(2:10:end) = 1;

dpar = m(:,:,:,circshift(ind,0,1));
dperp = m(:,:,:,circshift(ind,1,1));
theta = m(:,:,:,circshift(ind,2,1));
phi = m(:,:,:,circshift(ind,3,1));
d0 = m(:,:,:,circshift(ind,4,1));
rpar = m(:,:,:,circshift(ind,5,1));
rperp = m(:,:,:,circshift(ind,6,1));
r1 = m(:,:,:,circshift(ind,7,1));
r2 = m(:,:,:,circshift(ind,8,1));
w = m(:,:,:,circshift(ind,9,1));
