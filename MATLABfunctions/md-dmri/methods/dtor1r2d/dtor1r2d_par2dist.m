function dtor1r2d = dtor1r2d_par2dist(par,perp,theta,phi,d0,rpar,rperp,r1,r2,w)

n = numel(par);

if n>0
    dtor1r2d = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; r1'; r2'; w'];
    dtor1r2d = [n; dtor1r2d(:)];
else
    dtor1r2d = [];
end