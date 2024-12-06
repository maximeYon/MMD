function res = pov_tri2povinc(tri,povinc_path)

[parentdir,newdir] = fileparts(povinc_path);
if ~isfolder(povinc_path)
    mkdir(parentdir,newdir);
end

if ~isfield(tri,'n')
    tri.n = size(tri.verts,1);
end

if ~isfield(tri,'c')
    tri.c = ones(tri.n,3);
end

if ~isfield(tri,'norms')
    tri.norms = vertexNormal(triangulation(tri.tri, tri.verts),(1:tri.n)');
end

fnam = fullfile(povinc_path,'verts.txt');
fid = fopen(fnam, 'w');
N = size(tri.verts, 1);
fprintf(fid, '%8.0i,\n', N);
format = '%8.5f, %8.5f, %8.5f,\n';
for n = 1:N
    fprintf(fid,format,tri.verts(n,:));
end
fclose(fid);

fnam = fullfile(povinc_path,'tri.txt');
fid = fopen(fnam, 'w');
N = size(tri.tri, 1);
fprintf(fid, '%8.0i,\n', N);
format = '%8.0f, %8.0f, %8.0f,\n';
for n = 1:N
    fprintf(fid,format,tri.tri(n,:)-1);
end
fclose(fid);    

fnam = fullfile(povinc_path,'norms.txt');
fid = fopen(fnam, 'w');
N = size(tri.norms, 1);
fprintf(fid, '%8.0i,\n', N);
format = '%8.5f, %8.5f, %8.5f,\n';
for n = 1:N
    fprintf(fid,format,tri.norms(n,:));
end
fclose(fid);

fnam = fullfile(povinc_path,'c.txt');
fid = fopen(fnam, 'w');
N = size(tri.c, 1);
fprintf(fid, '%8.0i,\n', N);
format = '%8.5f, %8.5f, %8.5f,\n';
for n = 1:N
    fprintf(fid,format,tri.c(n,:));
end
fclose(fid);

if isfield(tri,'edges')
    fnam = fullfile(povinc_path,'edges.txt');
    fid = fopen(fnam, 'w');
    N = size(tri.edges, 1);
    fprintf(fid, '%8.0i,\n', N);
    format = '%8.0f, %8.0f, \n';
    for n = 1:N
        fprintf(fid,format,tri.edges(n,:)-1);
    end
    fclose(fid);
end


res = 1;
