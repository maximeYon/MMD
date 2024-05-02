function s = vasco16_1d_fit2data(m, xps)
% function s = vasco16_1d_fit2data(m, xps)


% convert to readable parameters
s0       = m(1);
v2       = m(2);
f_blood  = m(3);
D_blood  = m(4);
D_tissue = m(5);

% calculate tissue components and final fit 
e_tissue = exp(-xps.b * D_tissue);
e_blood  = exp(-xps.b * D_blood - xps.alpha2 * v2);

s = s0 * ( (1 - f_blood) * (e_tissue) + f_blood * (e_blood));

% allow different weights for different signal acqs
s(xps.s_ind == 1) = s(xps.s_ind == 1) * m(6);
