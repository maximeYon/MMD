function s = dtd_gamma_1d_fit2data(m, xps)
% function s = dtd_gamma_1d_fit2data(m, xps)
%
% s0        = m(1);
% d_iso     = m(2);
% mu2_iso   = m(3);
% mu2_aniso = m(4);
% 

% convert to readable parameters
s0        = m(1);
d_iso     = m(2);
mu2_iso   = m(3);
mu2_aniso = m(4);

mu2 = mu2_iso + mu2_aniso.*xps.b_delta.^2.*(xps.b_eta.^2+3)/3;

s = s0.*(1 + xps.b.*mu2./d_iso).^(-1*(d_iso.^2./mu2));        

