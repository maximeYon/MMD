function [n,par,perp,theta,phi,d0,rpar,rperp,w] = dtod_dist2par(dtod)

n = dtod(1);

if n>0
    m = numel(dtod(2:end))/n;
    dtod_array = reshape(dtod(2:end),[m n]);
    par = dtod_array(1,:)';
    perp = dtod_array(2,:)';
    theta = dtod_array(3,:)';
    phi = dtod_array(4,:)';
    d0 = dtod_array(5,:)';
    rpar = dtod_array(6,:)';
    rperp = dtod_array(7,:)';
    w = dtod_array(8,:)';
else
    par = [];
    perp = [];
    theta = [];
    phi = [];
    d0 = [];
    rpar = [];
    rperp = [];
    w = [];
end
