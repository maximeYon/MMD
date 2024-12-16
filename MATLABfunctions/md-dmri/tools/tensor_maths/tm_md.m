function md = tm_md(d2)
% function md = tm_md(d2)
% 
% Calculate the mean diffusivity of the second-order tensor 'd2'

md = tm_inner(d2, tm_3x3_to_1x6(tm_3x3_iso()));