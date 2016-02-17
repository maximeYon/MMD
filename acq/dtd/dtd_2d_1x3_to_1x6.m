function t = ts_2d_1x3_to_1x6(t)

xx = t(:,1);
yy = t(:,2);
xy = t(:,3);


% [xx yy c*xy]' * [xx yy c*xy] = 
%
% 1*xxxx 1*xxyy c *xxxy
% 1*yyxx 1*yyyy c *yyxy
% c*xyxx c*xyyy c2*xyxy
%
% symmetry gives
%
% xxxx xxyy xxxy
% xxyy yyyy yyxy
% xxxy yyxy xyxy
%

t = [ ...
    [xx.*xx yy.*yy] ...
    [xx.*yy] * sqrt(2) ...
    [xx.*xy yy.*xy] * sqrt(2) ...
    [xy.*xy] ];
