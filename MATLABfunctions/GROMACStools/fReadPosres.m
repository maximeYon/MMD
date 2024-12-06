function struct = fReadPosres(fname)

format = '%d%d%f%f%f';

fid = fopen(fname, 'r');
text = textscan(fid, format, 'headerlines', 1,'CommentStyle',';');
fclose(fid);

struct.atom = text{1};
struct.type = text{2};
struct.fx = text{3};
struct.fy = text{4};
struct.fz = text{5};
