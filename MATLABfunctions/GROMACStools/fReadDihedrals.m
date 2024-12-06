function struct = fReadDihedrals(fname)

format = '%d%d%d%d%d%s';

fid = fopen(fname, 'r');
text = textscan(fid, format, 'headerlines', 1,'CommentStyle',';');
fclose(fid);

struct.atom1 = text{1};
struct.atom2 = text{2};
struct.atom3 = text{3};
struct.atom4 = text{4};
struct.function = text{5};
struct.gromostype = text{6};
