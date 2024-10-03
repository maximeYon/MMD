function gro_out = fSubsampleGro(gro_in,index)

gro_out = gro_in;
gro_out.Natom = length(index);
gro_out.molnr = gro_in.molnr(index);
gro_out.molecule = gro_in.molecule(index);
gro_out.atom = gro_in.atom(index);
gro_out.atomnr = gro_in.atomnr(index);
gro_out.x = gro_in.x(index);
gro_out.y = gro_in.y(index);
gro_out.z = gro_in.z(index);