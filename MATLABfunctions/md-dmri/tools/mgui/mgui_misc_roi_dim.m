function [c_dim, w_out, n, w] = mgui_misc_roi_dim(I_roi, n)
% function c_dim = mgui_misc_roi_dim(I_roi)

if (nargin < 2)
    n = sum(I_roi(:) > 0);
end

w = [];

if (n == 0)
    c_dim = [];
    w_out = [];
    n = [];
    w = [];
    return;
end

w = zeros(1,3);
w(1) = sum(sum(sum(I_roi,2),3) > 0);
w(2) = sum(sum(sum(I_roi,1),3) > 0);
w(3) = sum(sum(sum(I_roi,1),2) > 0);

[w_out, c_dim] = min(w);



