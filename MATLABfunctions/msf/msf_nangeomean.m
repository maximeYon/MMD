function a = msf_nangeomean(a,d)
% function a = msf_nangeomean(a,d)

i = isnan(a);
a(i) = 1;
a = sum(log(a),d) ./ sum(~i,d);
a = exp(a);
a(a == Inf) = 0;