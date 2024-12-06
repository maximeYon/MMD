function gro1 = fNatomGro(gro,atomindex)

gro1 = gro;
gro1.Natom = length(atomindex);
gro1.molnr = gro.molnr(atomindex);
gro1.molecule = gro.molecule(atomindex);
gro1.atom = gro.atom(atomindex);
gro1.atomnr = gro.atomnr(atomindex);
gro1.x = gro.x(atomindex);
gro1.y = gro.y(atomindex);
gro1.z = gro.z(atomindex);