function pairs = fReadPairs(fname)

format = '%d%d%d';

fid = fopen(fname, 'r');
text = textscan(fid, format, 'headerlines', 1,'CommentStyle',';');
fclose(fid);

pairs.atom1 = text{1};
pairs.atom2 = text{2};
pairs.function = text{3};

