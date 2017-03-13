function [E_bulk, E_shear, E2_iso] = tm_1x21_iso()
% function [E_bulk, E_shear, E2_iso] = tm_1x21_iso()
% 
% Obtain the isotropic fourth-order tensors in the 1x21 Voigt-like format


[E_bulk, E_shear, E2_iso] = tm_6x6_iso;
E_bulk  = tm_6x6_to_1x21(E_bulk);
E_shear = tm_6x6_to_1x21(E_shear);
E2_iso  = tm_6x6_to_1x21(E2_iso);
