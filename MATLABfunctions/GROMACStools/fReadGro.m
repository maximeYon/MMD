function gro = fReadGro(fname,NoVelocity)

if nargin == 2
    %not used anymore
    if NoVelocity == 1
        CharPerLine = 45;
    else
        CharPerLine = 69;
    end
end


fid = fopen(fname, 'r');

for n=1:3
    templine = fgets(fid);
end
CharPerLine = length(templine);
frewind(fid)

gro.header = fgetl(fid);
text.Natom = textscan(fid, '%5d', 1,'headerlines',0);
gro.Natom = text.Natom{1};

text = fscanf(fid,'%c',double(gro.Natom*CharPerLine));
text = reshape(text,CharPerLine,gro.Natom)';
gro.molnr = str2num(text(:,1+(1:5)));
gro.molecule = cellstr(text(:,6+(1:5)));
gro.atom = deblank(strjust(cellstr(text(:,11+(1:5))),'left'));
gro.atomnr = str2num(text(:,16+(1:5)));
gro.x = str2num(text(:,21+(1:8)));
gro.y = str2num(text(:,29+(1:8)));
gro.z = str2num(text(:,37+(1:8)));

format = '%10.5f%10.5f%10.5f';
box = textscan(fid, format, 1);
fclose(fid);

gro.boxxx = box{1};
gro.boxyy = box{2};
gro.boxzz = box{3};