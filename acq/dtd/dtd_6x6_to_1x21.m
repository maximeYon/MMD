function t = ts_6x6_to_1x21(t)


% xxxx xxyy xxzz c1*xxxy c1*xxxz c1*xxyz
% xxyy xxxx yyzz c1*yyyx c1*yyxz c1*yyyz
% xxzz yyzz zzzz c1*zzxy c1*zzzx c1*zzzy
%  -    -    -   c2*xyxy c2*xyxz c2*xyyz
%  -    -    -   c2*xyxz c2*xzxz c2*xzyz
%  -    -    -   c2*xyyz c2*xzyz c2*yzyz


t = [...
    [t(1,1) t(2,2) t(3,3)] * sqrt(1) ...
    [t(1,2) t(1,3) t(2,3)] * sqrt(2) ...
    [t(1,6) t(5,2) t(4,3)] * sqrt(2) ... % xxyz yyxz zzxy
    [t(1,4) t(1,5) t(2,4)] * sqrt(2) ...
    [t(2,6) t(3,5) t(3,6)] * sqrt(2) ...
    [t(4,4) t(5,5) t(6,6)] * sqrt(1) ...
    [t(4,5) t(4,6) t(5,6)] * sqrt(2)];

% 
% xxxx = t(1); % all same
% yyyy = t(2);
% zzzz = t(3);
% 
% c = sqrt(4);
% xxyy = t(4) / c; % two same
% xxzz = t(5) / c;
% yyzz = t(6) / c;
% 
% c = sqrt( 2 );
% xxyz = t(7) / c; % three species
% yyxz = t(8) / c;
% zzxy = t(9) / c;
% 
% c = sqrt( 2 );
% xxxy = t(10) / c; % two species, three same
% xxxz = t(11) / c;
% yyxy = t(12) / c;
% yyyz = t(13) / c;
% zzxz = t(14) / c;
% zzyz = t(15) / c;
% 
% c = 1;
% xyxy = t(16) / c;
% xzxz = t(17) / c;
% yzyz = t(18) / c;
% 
% c = sqrt(2);
% xyxz = t(19) / c;
% xyyz = t(20) / c;
% xzyz = t(21) / c;
% 
% 
% A = [...
%     xxxx xxyy xxzz; 
%     xxyy yyyy yyzz;
%     xxzz yyzz zzzz]; % 6 unique
% 
% B = [...
%     xxxy xxxz xxyz;
%     yyxy yyxz yyyz;
%     zzxy zzxz zzyz]; % 9 unique
% 
% C = [...
%     xyxy xyxz xyyz; 
%     xyxz xzxz xzyz; 
%     xyyz xzyz yzyz]; 
%     
% % this is fishy, but probably correct given the assumptions on cross terms
% % w = sqrt(ones(3,3)*2 + eye(3)*2);
% % w = ones(3,3);
% 
% t = [A B; B' C];
%     
