function [E_bulk, E_shear, E2_iso] = dtd_1x21_iso()
% function [E_bulk, E_shear, E2_iso] = dtd_1x21_iso()
% 
% Obtain the isotropic fourth-order tensors in the 1x21 Voigt-like format


[E_bulk, E_shear, E2_iso] = ts_6x6_iso;
E_bulk = ts_6x6_to_1x21(E_bulk);
E_shear = ts_6x6_to_1x21(E_shear);
