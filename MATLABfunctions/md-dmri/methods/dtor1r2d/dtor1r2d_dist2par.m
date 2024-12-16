function [n,par,perp,theta,phi,d0,rpar,rperp,r1,r2,w] = dtor1r2d_dist2par(dtor1r2d)

n = dtor1r2d(1);

if n>0
    m = numel(dtor1r2d(2:end))/n;
    dtor1r2d_array = reshape(dtor1r2d(2:end),[m n]);
    par = dtor1r2d_array(1,:)';
    perp = dtor1r2d_array(2,:)';
    theta = dtor1r2d_array(3,:)';
    phi = dtor1r2d_array(4,:)';
    d0 = dtor1r2d_array(5,:)';
    rpar = dtor1r2d_array(6,:)';
    rperp = dtor1r2d_array(7,:)';
    r1 = dtor1r2d_array(8,:)';
    r2 = dtor1r2d_array(9,:)';
    w = dtor1r2d_array(10,:)';
else
    par = [];
    perp = [];
    theta = [];
    phi = [];
    d0 = [];
    rpar = [];
    rperp = [];
    r1 = [];
    r2 = [];
    w = [];
end
