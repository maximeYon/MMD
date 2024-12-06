function phi = fHydrogenOrbital(CartCoord,Para,QuantNum);

x = CartCoord.x;
y = CartCoord.y;
z = CartCoord.z;

a0 = Para.a0;
Z = Para.Z;

Principal = QuantNum.n;
OrbitalAngularMomentum = QuantNum.l;
Magnetic = QuantNum.m;

r = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z./r);
phi = angle(x + i*y);

if strcmp(Principal,'1') == 1
    R = 2*(Z/a0)^1.5*exp(-Z*r/a0);
elseif strcmp(Principal,'2') == 1
    if strcmp(OrbitalAngularMomentum,'s') == 1
        R = 1/2/sqrt(2)*(Z/a0)^1.5*(2-Z*r/a0).*exp(-Z*r/2/a0);
    elseif strcmp(OrbitalAngularMomentum,'p') == 1
        R = 1/2/sqrt(6)*(Z/a0)^1.5*(Z*r/a0).*exp(-Z*r/2/a0);
    end
elseif strcmp(Principal,'3') == 1
    if strcmp(OrbitalAngularMomentum,'s') == 1
        R = 2/9/sqrt(3)*(Z/a0)^1.5*(3-2*Z*r/a0+2/9*(Z*r/a0).^2).*exp(-Z*r/3/a0);
    elseif strcmp(OrbitalAngularMomentum,'p') == 1
        R = 2/9/sqrt(6)*(Z/a0)^1.5*(2-Z*r/3/a0).*exp(-Z*r/3/a0);
    elseif strcmp(OrbitalAngularMomentum,'d') == 1
        R = 4/81/sqrt(30)*(Z/a0)^1.5*(Z*r/a0).^2.*exp(-Z*r/3/a0);
    end
end

if strcmp(Magnetic,'0') == 1
    Y = sqrt(1/4/pi);
elseif strcmp(Magnetic,'x') == 1
    Y = sqrt(3/4/pi)*sin(theta).*cos(phi);
elseif strcmp(Magnetic,'y') == 1
    Y = sqrt(3/4/pi)*sin(theta).*sin(phi);
elseif strcmp(Magnetic,'z') == 1
    Y = sqrt(3/4/pi)*cos(theta);
elseif strcmp(Magnetic,'xy') == 1
    Y = sqrt(15/16/pi)*sin(theta).^2.*sin(2*phi);
elseif strcmp(Magnetic,'yz') == 1
    Y = sqrt(15/4/pi)*cos(theta).*sin(theta).*sin(phi);
elseif strcmp(Magnetic,'zx') == 1
    Y = sqrt(15/4/pi)*cos(theta).*sin(theta).*cos(phi);
elseif strcmp(Magnetic,'x2-y2') == 1
    Y = sqrt(15/16/pi)*sin(theta).^2.*cos(2*phi);
elseif strcmp(Magnetic,'z2') == 1
    Y = sqrt(5/16/pi)*(3*cos(theta).^2-1);
end

phi = R.*Y;