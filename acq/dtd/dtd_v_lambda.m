function v_lambda = ts_v_lambda(d2)

[~,E_shear] = ts_6x6_iso;

v_lambda = ts_inner(ts_1x6_to_1x21(d2), ts_6x6_to_1x21(E_shear));