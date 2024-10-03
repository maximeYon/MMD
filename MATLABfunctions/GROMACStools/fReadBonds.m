function bonds = fReadBonds(fname)

format = '%d%d%d%s';

fid = fopen(fname, 'r');
text = textscan(fid, format, 'headerlines', 1,'CommentStyle',';');
fclose(fid);

bonds.atom1 = text{1};
bonds.atom2 = text{2};
bonds.function = text{3};
bonds.gromostype = text{4};

