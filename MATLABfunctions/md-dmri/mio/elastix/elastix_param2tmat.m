function Trafo = elastix_param2tmat(Par)
% function Trafo = elastix_param2tmat(Par)
%
% From Alexander Leemans 2011
%
% Ex: Trafo = E_DTI_ConvertP2TMat([0 0 0; 1 2 3; 1 1 1; 0 0 0]);
% [rotx roty rotz; transx transy transz; scalex scaley scalez; skewx skewy
% skewz];
%
% roti in radians

if (numel(Par) ~= 12)
    error('Requires 12 parameters');
elseif (size(Par,1) == 1)
    Par = reshape(Par, [3 4])';
elseif (size(Par,1) ~= 4)
    error('Size must be 4x3 or 1x12');
end

tx = Par(2,1);
ty = Par(2,2);
tz = Par(2,3);

rx = Par(1,1);
ry = Par(1,2);
rz = Par(1,3);

sx = 1/Par(3,1);
sy = 1/Par(3,2);
sz = 1/Par(3,3);

gx = Par(4,1);
gy = Par(4,2);
gz = Par(4,3);

T = [1 0 0 tx;
    0 1 0 ty;
    0 0 1 tz;
    0 0 0 1];

Rx = [1 0 0 0;
    0 cos(rx) sin(rx) 0;
    0 -sin(rx) cos(rx) 0;
    0 0 0 1];

Ry = [cos(ry) 0 -sin(ry) 0;
    0 1 0 0;
    sin(ry) 0 cos(ry) 0;
    0 0 0 1];

Rz = [cos(rz) sin(rz) 0 0;
    -sin(rz) cos(rz) 0 0;
    0 0 1 0;
    0 0 0 1];

R = Rx*Ry*Rz;

Gx = [1 0 0 0;
    0 1 0 0;
    -gx 0 1 0;
    0 0 0 1];

Gy = [1 -gy 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

Gz = [1 0 0 0;
    0 1 -gz 0;
    0 0 1 0;
    0 0 0 1];

G = Gx*Gy*Gz;

S = [sx 0 0 0;
    0 sy 0 0;
    0 0 sz 0;
    0 0 0 1];


Trafo = T*inv(R*G*S);


