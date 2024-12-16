function dtod = dtod_par2dist(par,perp,theta,phi,d0,rpar,rperp,w)

n = numel(par);

if n>0
    dtod = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; w'];
    dtod = [n; dtod(:)];
else
    dtod = [];
end