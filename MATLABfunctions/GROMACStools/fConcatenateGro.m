function gro_out = fConcatenateGro(gro1,gro2)

gro_out = gro1;

gro_out.Natom = numel(gro1.x) + numel(gro2.x);
gro_out.molnr = cat(1,gro1.molnr,max(gro1.molnr)+gro2.molnr);
gro_out.molecule = cat(1,gro1.molecule,gro2.molecule);
gro_out.atom = cat(1,gro1.atom,gro2.atom);
gro_out.atomnr = cat(1,gro1.atomnr,max(gro1.atomnr)+gro2.atomnr);
gro_out.x = cat(1,gro1.x,gro2.x);
gro_out.y = cat(1,gro1.y,gro2.y);
gro_out.z = cat(1,gro1.z,gro2.z);