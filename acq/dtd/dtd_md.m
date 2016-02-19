function md = dtd_md(d2)
% function md = dtd_md(d2)
% 
% Calculate the mean diffusivity of the second-order tensor 'd2'

md = dtd_inner(d2, dtd_3x3_to_1x6(dtd_3x3_iso()));