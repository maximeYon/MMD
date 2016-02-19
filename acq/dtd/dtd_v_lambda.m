function v_lambda = dtd_v_lambda(d2)
% function v_lambda = dtd_v_lambda(d2)
%
% Calculate the eigenvalue variance of the second-order tensor in 'd2'

[~,E_shear] = dtd_6x6_iso;

v_lambda = dtd_inner(dtd_1x6_to_1x21(d2), dtd_6x6_to_1x21(E_shear));