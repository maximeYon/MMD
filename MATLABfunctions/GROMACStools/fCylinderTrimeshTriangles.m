function triangles = fCylinderTrimeshTriangles(Nphi,Nz)

triangles = [];
for nz = 1:(Nz-1)
    triangles_temp = [(1:Nphi)' circshift((1:Nphi)',-1,1) (1:Nphi)'+Nphi] + (nz-1)*Nphi;
    triangles = cat(1,triangles,triangles_temp);
end
for nz = 1:(Nz-1)
    triangles_temp = [(1:Nphi)' (1:Nphi)'+Nphi circshift((1:Nphi)'+Nphi,1,1) ] + (nz-1)*Nphi;
    triangles = cat(1,triangles,triangles_temp);
end

