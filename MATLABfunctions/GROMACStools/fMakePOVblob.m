function res = fMakePOVblob(fnam,gro,molnr,chainatoms)

fid = fopen([fnam '.inc'],'w');

for ncount = 1:length(molnr)
    n = molnr(ncount);
    index = [];
    for m = 1:length(chainatoms)
        index = [index; find(all([strcmp(gro.atom, chainatoms{m}) gro.molnr==n],2))];
    end
    
    indexmid = round(length(index)/2);
    
    indexmove = find((gro.x(index)-gro.x(index(indexmid))) > .5*gro.boxxx);
    gro.x(index(indexmove)) = gro.x(index(indexmove)) - gro.boxxx;
    indexmove = find((gro.x(index)-gro.x(index(indexmid))) < -.5*gro.boxxx);
    gro.x(index(indexmove)) = gro.x(index(indexmove)) + gro.boxxx;
    indexmove = find((gro.y(index)-gro.y(index(indexmid))) > .5*gro.boxyy);
    gro.y(index(indexmove)) = gro.y(index(indexmove)) - gro.boxyy;
    indexmove = find((gro.y(index)-gro.y(index(indexmid))) < -.5*gro.boxyy);
    gro.y(index(indexmove)) = gro.y(index(indexmove)) + gro.boxyy;
    indexmove = find((gro.z(index)-gro.z(index(indexmid))) > .5*gro.boxzz);
    gro.z(index(indexmove)) = gro.z(index(indexmove)) - gro.boxzz;
    indexmove = find((gro.z(index)-gro.z(index(indexmid))) < -.5*gro.boxzz);
    gro.z(index(indexmove)) = gro.z(index(indexmove)) + gro.boxzz;
    
    tempstr = ['blob{']; 
    fprintf(fid,'%s\n',tempstr);
    for m = 1:length(index)
        tempstr = ['sphere{<' num2str(gro.x(index(m))) ','...
            num2str(gro.z(index(m))) ',' ...
            num2str(gro.y(index(m))) '>' ',blobrad,blobstrength}']; 
        fprintf(fid,'%s\n',tempstr);
    end
    tempstr = ['threshold blobthresh}']; 
    fprintf(fid,'%s\n\n',tempstr);
end
fclose(fid);

res = 1;