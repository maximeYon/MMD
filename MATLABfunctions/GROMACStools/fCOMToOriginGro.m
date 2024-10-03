function gro_out = fCOMToOriginGro(gro_in)

gro = gro_in;

gro.x = gro.x - mean(gro.x);
gro.y = gro.y - mean(gro.y);
gro.z = gro.z - mean(gro.z);

gro_out = gro;