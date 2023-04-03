function angles = fReadAngles(fname)

format = '%d%d%d%d%s';

fid = fopen(fname, 'r');
text = textscan(fid, format, 'headerlines', 1,'CommentStyle',';');
fclose(fid);

angles.atom1 = text{1};
angles.atom2 = text{2};
angles.atom3 = text{3};
angles.function = text{4};
angles.gromostype = text{5};
