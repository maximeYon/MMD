function t = dtd_1x6_to_1x21(n2)
% function t = dtd_1x6_to_1x21(n2)
%
% Convert a second-order tensor on Voigt notation (1x6) to a 
% fourth-order tensor in Voight-like notation (1x21). This correspons
% to taking the outer product of 'n2'

xx = n2(:,1);
yy = n2(:,2);
zz = n2(:,3);
xy = n2(:,4);
xz = n2(:,5);
yz = n2(:,6);

% % Kurtosis tensor direction matrix
% if (~do_limit)
%     tmp = repmat([1 1 1 6 6 6 4 4 4 4 4 4 12 12 12],size(n,1),1);
%     n_KT = [ ...
%         tmp .* ...
%         n(:,[1 2 3 1 2 1 1 2 1 1 2 1 1 1 1]) .* ...
%         n(:,[1 2 3 1 2 1 1 2 3 2 3 1 2 1 2]) .* ...
%         n(:,[1 2 3 2 3 3 1 2 3 2 3 1 3 2 2]) .* ...
%         n(:,[1 2 3 2 3 3 2 3 3 2 3 3 3 3 3])];
% else


% [xx yy zz c*xy c*xz c*yz]' * [xx yy zz c*xy c*xz c*yz]
%
% xxxx xxyy xxzz c1*xxxy c1*xxxz c1*xxyz
% xxyy xxxx yyzz c1*yyyx c1*yyxz c1*yyyz
% xxzz yyzz zzzz c1*zzxy c1*zzzx c1*zzzy
%  -    -    -   c2*xyxy c2*xyxz c2*xyyz
%  -    -    -   c2*xyxz c2*xzxz c2*xzyz
%  -    -    -   c2*xyyz c2*xzyz c2*yzyz
% 



t = [...
    [xx.*xx yy.*yy zz.*zz] * sqrt(1) ...     % 1 = 1 *
    [xx.*yy xx.*zz yy.*zz] * sqrt(2) ...     % 1 = 2
    [xx.*yz yy.*xz zz.*xy] * sqrt(2) ...     % 2 = 4
    [xx.*xy xx.*xz]        * sqrt(2) ...     % 2 = 4
    [yy.*xy yy.*yz]        * sqrt(2) ...     % 2 = 4
    [zz.*xz zz.*yz]        * sqrt(2) ...     % 2 = 4
    [xy.*xy xz.*xz yz.*yz] * sqrt(1) ...     % 4 = 4
    [xy.*xz xy.*yz xz.*yz] * sqrt(2) ];      % 4 = 8


    