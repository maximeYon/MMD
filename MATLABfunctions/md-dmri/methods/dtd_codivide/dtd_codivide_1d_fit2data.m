function s = dtd_codivide_1d_fit2data(m, xps)
% function s = dtd_codivide_1d_fit2data(m, xps)
%
% s0     = m(1);
% v_at   = m(2);
% v_fw   = m(3);
% md     = m(4);
% md_fw  = m(5);
% 

s0     = m(1);
v_at   = m(2);
v_fw   = m(3);
md     = m(4);
md_fw  = m(5);
sw     = m(6:end);


v_it   = 1 - v_at - v_fw;

ad_at = md * 3;
rd_at = 0;

md_it = md;


md_at      = (ad_at + 2 * rd_at) / 3;
d_delta_at = (ad_at - rd_at) / (3 * md_at);

a_at  = 3*xps.b.*md_at.*xps.b_delta.*d_delta_at + eps;

s_at = exp(-xps.b.*md_at).*exp(a_at/3).*...
    sqrt(pi)/2.*real(gammainc(a_at,1/2)./sqrt(a_at));


s_it = exp(-xps.b * md_it);


s_fw = exp(-xps.b * md_fw);



s = s0 * ( ...
    v_fw * s_fw + ...
    v_at * s_at + ...
    v_it * s_it); 

for c = 1:numel(sw)
    ind = (xps.s_ind == (c + 1));
    s(ind) = s(ind) * sw(c);
end
    
