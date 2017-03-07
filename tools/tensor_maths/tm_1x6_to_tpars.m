function t = tm_1x6_to_tpars(t_1x6)
% function t = tm_1x6_to_tpars(t_1x6)

if (size(t_1x6,2) ~= 6), error('wrong dimension of tensor'); end


t = tm_3x3_to_tpars(tm_1x6_to_3x3(t_1x6));


