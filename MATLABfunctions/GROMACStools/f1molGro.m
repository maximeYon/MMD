function gro1 = f1molGro(gro,molnr)

index = find(gro.molnr == molnr);

gro1 = gro;
gro1.Natom = length(index);
gro1.molnr = gro.molnr(index);
gro1.molecule = gro.molecule(index);
gro1.atom = gro.atom(index);
gro1.atomnr = gro.atomnr(index);
gro1.x = gro.x(index);
gro1.y = gro.y(index);
gro1.z = gro.z(index);