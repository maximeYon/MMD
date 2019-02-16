function a = msf_nanmean(a,d)
% function a = msf_nanmean(a,d)

if nargin < 2
    d = 1;
end

i = isnan(a);
a(i) = 0;
a = sum(a,d) ./ sum(~i,d);
a(a == Inf) = 0;