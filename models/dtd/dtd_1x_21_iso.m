function [E_bulk, E_shear, E2_iso] = ts_1x21_iso()


[E_bulk, E_shear, E2_iso] = ts_6x6_iso;
E_bulk = ts_6x6_to_1x21(E_bulk);
E_shear = ts_6x6_to_1x21(E_shear);
