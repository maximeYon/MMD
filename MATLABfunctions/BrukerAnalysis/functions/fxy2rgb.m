function c = fxy2rgb(X,Y,Size)

mid.x = Size(1); mid.y = Size(2);
width.x = Size(3); width.y = Size(4);


Z = 2*(X-mid.x)/width.x + i*2*(Y-mid.y)/width.y;
R = abs(Z);
phi = angle(Z)/2/pi;
phi = phi - 1/6;

index = find(phi<0);
phi(index) = phi(index)+1;

ctrans.r = zeros(size(phi));
ctrans.g = zeros(size(phi));
ctrans.b = zeros(size(phi));


index = find(phi>=0 & phi<1/6);
ctrans.r(index) = 1;
ctrans.g(index) = phi(index)*6;

index = find(phi>=1/6 & phi<1/3);
ctrans.g(index) = 1;
ctrans.r(index) = 1-(phi(index)-1/6)*6;

index = find(phi>=1/3 & phi<1/2);
ctrans.g(index) = 1;
ctrans.b(index) = (phi(index)-1/3)*6;

index = find(phi>=1/2 & phi<2/3);
ctrans.b(index) = 1;
ctrans.g(index) = 1-(phi(index)-1/2)*6;

index = find(phi>=2/3 & phi<5/6);
ctrans.b(index) = 1;
ctrans.r(index) = (phi(index)-2/3)*6;

index = find(phi>=5/6 & phi<1);
ctrans.r(index) = 1;
ctrans.b(index) = 1-(phi(index)-5/6)*6;

R(find(R>1)) = 1;

c.r = R.*ctrans.r + (1-R);
c.g = R.*ctrans.g + (1-R);
c.b = R.*ctrans.b + (1-R);
    