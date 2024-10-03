function res = fWriteGro(fname,gro)

fid = fopen(fname,'w');

fprintf(fid,'%s\n',gro.header);
fprintf(fid,'%5d\n',gro.Natom);

format = '%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n';

for natom = 1:gro.Natom;
    fprintf(fid,format,gro.molnr(natom),gro.molecule{natom},gro.atom{natom},...
        gro.atomnr(natom),gro.x(natom),gro.y(natom),gro.z(natom));
end

format = '%10.5f%10.5f%10.5f\n';
fprintf(fid,format,gro.boxxx,gro.boxyy,gro.boxzz);

fclose(fid);

res = 1;