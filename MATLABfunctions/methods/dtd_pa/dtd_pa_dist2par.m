function [n,par,perp,w] = dtd_pa_dist2par(dtd_pa)

n = dtd_pa(1);

if n>0
    m = numel(dtd_pa(2:end))/n;
    dtd_array = reshape(dtd_pa(2:end),[m n]);
    par = dtd_array(1,:)';
    perp = dtd_array(2,:)';
    w = dtd_array(3,:)';
else
    par = [];
    perp = [];
    w = [];
end
