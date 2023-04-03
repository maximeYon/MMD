function GyrationTensor = fGyrationTensorGro(gro)

gro = fCOMToOriginGro(gro);

Natoms = numel(gro.x);
xx = 1/Natoms*sum(gro.x.*gro.x);
yy = 1/Natoms*sum(gro.y.*gro.y);
zz = 1/Natoms*sum(gro.z.*gro.z);
xy = 1/Natoms*sum(gro.x.*gro.y);
xz = 1/Natoms*sum(gro.x.*gro.z);
yz = 1/Natoms*sum(gro.y.*gro.z);

GyrationTensor = tm_1x6_to_tpars([xx yy zz sqrt(2)*[xy xz yz]]);
