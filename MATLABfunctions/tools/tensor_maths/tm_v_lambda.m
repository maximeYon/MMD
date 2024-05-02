function v_lambda = tm_v_lambda(d2)
% function v_lambda = tm_v_lambda(d2)
%
% Calculate the eigenvalue variance of the second-order tensor in 'd2'

[~,E_shear] = tm_6x6_iso;

v_lambda = tm_inner(tm_1x6_to_1x21(d2), tm_6x6_to_1x21(E_shear));