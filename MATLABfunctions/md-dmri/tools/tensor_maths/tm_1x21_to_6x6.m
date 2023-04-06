function t = tm_1x21_to_6x6(t)
% function t = tm_1x21_to_6x6(t)
% 
% Convert fourth-order tensor in 1x21 format to 6x6 format

xxxx = t(1); % all same
yyyy = t(2);
zzzz = t(3);

c = sqrt( 2 );
xxyy = t(4) / c; 
xxzz = t(5) / c;
yyzz = t(6) / c;

c = sqrt( 2 );
xxyz = t(7) / c; 
yyxz = t(8) / c;
zzxy = t(9) / c;

c = sqrt( 2 );
xxxy = t(10) / c;
xxxz = t(11) / c;
yyxy = t(12) / c;
yyyz = t(13) / c;
zzxz = t(14) / c;
zzyz = t(15) / c;

c = 1;
xyxy = t(16) / c;
xzxz = t(17) / c;
yzyz = t(18) / c;

c = sqrt(2);
xyxz = t(19) / c;
xyyz = t(20) / c;
xzyz = t(21) / c;


A = [...
    xxxx xxyy xxzz; 
    xxyy yyyy yyzz;
    xxzz yyzz zzzz]; % 6 unique

B = [...
    xxxy xxxz xxyz;
    yyxy yyxz yyyz;
    zzxy zzxz zzyz]; % 9 unique

C = [...
    xyxy xyxz xyyz; 
    xyxz xzxz xzyz; 
    xyyz xzyz yzyz]; 
    
% this is fishy, but probably correct given the assumptions on cross terms
% w = sqrt(ones(3,3)*2 + eye(3)*2);
% w = ones(3,3);

t = [A B; B' C];
    
