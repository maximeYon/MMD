function res = fSpheres2POVinc(fnam,Spheres)

fid = fopen(fnam,'w');

Nnodes = length(Spheres.x);
for nnode = 1:Nnodes
    tempstr = ['sphere{<' num2str(Spheres.x(nnode)) ','...
        num2str(Spheres.z(nnode)) ',' ...
        num2str(Spheres.y(nnode)) '>' ',rad}']; 
    fprintf(fid,'%s\n',tempstr);
end

fclose(fid);

res = 1;