function [E_bulk, E_shear, E2_iso] = ts_6x6_iso()

E2_iso = ts_3x3_to_1x6(eye(3)/3);

E_bulk = ts_outer(E2_iso, E2_iso);
E_shear = eye(6)/3 - E_bulk;