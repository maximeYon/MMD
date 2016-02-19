function [E_bulk, E_shear, E2_iso] = dtd_6x6_iso()
% function [E_bulk, E_shear, E2_iso] = dtd_6x6_iso()
%
% Get the three different isotropic fourth-order tensors

E2_iso = dtd_3x3_to_1x6(eye(3)/3);

E_bulk = dtd_outer(E2_iso, E2_iso);
E_shear = eye(6)/3 - E_bulk;