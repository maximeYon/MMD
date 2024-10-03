function res = fWriteNdx(fname,gro)

fid = fopen(fname,'w');

fprintf(fid,'%s\n','[ System ]');

format = '%8d\n';
for count = 1:gro.Natom
    fprintf(fid,format,count);
end

fclose(fid);

res = 1;
