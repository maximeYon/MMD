function fa = ts_fa(a,b)

if (nargin == 1)
    
    d2 = a;
    md = ts_md(d2);
    vl = ts_v_lambda(d2);
    
elseif (nargin == 2)
    
    md = a;
    vl = b;
    
end


fa = sqrt( (3/2) * (1 + md.^2./vl).^(-1) );