function fa = dtd_fa(a,b)
% function fa = dtd_fa(a,b)
%
% With one input argumnet: calculate the FA of the second order tensor 'a'
% With two arguments: assume 'a' to be MD and 'b' to be the variance in
%           eigenvectors, and calculate FA from that


if (nargin == 1)
    
    d2 = a;
    md = ts_md(d2);
    vl = ts_v_lambda(d2);
    
elseif (nargin == 2)
    
    md = a;
    vl = b;
    
end


fa = sqrt( (3/2) * (1 + md.^2./vl).^(-1) );