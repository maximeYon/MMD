function c = fmxyz2color(mxall,myall,mzall)

[Nt,Nz] = size(mxall);

for n = 1:3
    clong(:,:,n) = mzall'/2 + .5;
end

c = clong;

ctrans = zeros(size(clong));

phi = reshape(angle(mxall'+i*(myall+1e-10)')/2/pi,Nt*Nz,1);
pos = find(phi<0);
phi(pos) = phi(pos)+1;

ctransr = zeros(size(phi));
ctransg = zeros(size(phi));
ctransb = zeros(size(phi));


pos = find(phi>=0 & phi<1/6);
ctransr(pos) = 1;
ctransg(pos) = phi(pos)*6;

pos = find(phi>=1/6 & phi<1/3);
ctransg(pos) = 1;
ctransr(pos) = 1-(phi(pos)-1/6)*6;

pos = find(phi>=1/3 & phi<1/2);
ctransg(pos) = 1;
ctransb(pos) = (phi(pos)-1/3)*6;

pos = find(phi>=1/2 & phi<2/3);
ctransb(pos) = 1;
ctransg(pos) = 1-(phi(pos)-1/2)*6;

pos = find(phi>=2/3 & phi<5/6);
ctransb(pos) = 1;
ctransr(pos) = (phi(pos)-2/3)*6;

pos = find(phi>=5/6 & phi<1);
ctransr(pos) = 1;
ctransb(pos) = 1-(phi(pos)-5/6)*6;


ctrans(:,:,1) = reshape(ctransr,Nz,Nt);
ctrans(:,:,2) = reshape(ctransg,Nz,Nt);
ctrans(:,:,3) = reshape(ctransb,Nz,Nt);

mtrans = abs(mxall'+i*myall');

for n = 1:3
    c(:,:,n) = mtrans.*ctrans(:,:,n) + (1-mtrans).*clong(:,:,n);
end