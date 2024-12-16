function res = bruker_mkshapefile(out_fn,g,title_str)
% function res = bruker_mkshapefile(out_fn,g,title_str)

if nargin < 3
    title_str = 'Gradient shape generated by Matlab';
end

res = -1;

[out_path,out_name,out_ext] = fileparts(out_fn);
if ~isdir(out_path)
    mkdir(out_path)
end

fid = fopen(out_fn,'w');

text.header = {
{['##TITLE= ' title_str]};
{'##JCAMP-DX= 5.00 Bruker JCAMP library'}
{'##DATA TYPE= Shape Data'}
{'##ORIGIN= Bruker Analytik GmbH'}
{'##OWNER= <nmrsu>'}
{'##DATE= xx'}
{'##TIME= xx'}
{'##MINX= 0'}
{'##MAXX= 1'}
{'##MINY= 0'}
{'##MAXY= 1'}
{'##$SHAPE_EXMODE= gradient'}
%{'##$SHAPE_EXMODE= Gradient'}
{'##$SHAPE_TOTROT= 0'}
{'##$SHAPE_BWFAC= 0'}
{'##$SHAPE_INTEGFAC= 0'}
{'##$SHAPE_MODE= 0'}};

[nlines, ~] = size(text.header);

for nline = 1:nlines
    fprintf(fid,'%s\n',text.header{nline}{1});
end

fprintf(fid,'%s\n',['##NPOINTS=' num2str(numel(g))]);
fprintf(fid,'%s\n',['##XYDATA= (X++(Y..Y))']);

for ng = 1:numel(g)
    fprintf(fid,'%f\n',g(ng));
end

fprintf(fid,'%s\n',['##END']);
fclose(fid);

res = 1;