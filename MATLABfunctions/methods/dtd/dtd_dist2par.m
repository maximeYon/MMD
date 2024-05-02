function [n,par,perp,theta,phi,w] = dtd_dist2par(dtd)

n = dtd(1);

if n>0
    m = numel(dtd(2:end))/n;
    dtd_array = reshape(dtd(2:end),[m n]);
    par = dtd_array(1,:)';
    perp = dtd_array(2,:)';
    theta = dtd_array(3,:)';
    phi = dtd_array(4,:)';
    w = dtd_array(5,:)';
else
    par = [];
    perp = [];
    theta = [];
    phi = [];
    w = [];
end
