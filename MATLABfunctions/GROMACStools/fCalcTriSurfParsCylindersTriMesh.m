function tri = fCalcTriSurfParsCylindersTriMesh(gro,cylpars)

Nphi = cylpars(1);
Nz = cylpars(2);

% gro = fPutMoleculesInBoxGro(gro);
% gro = fAlignGyrationTensorGro(gro);
% gro = fCOMToOriginGro(gro);        
% %gro.y = gro.y - min(gro.y);

triangles = fCylinderTrimeshTriangles(Nphi,Nz);

molnr = unique(gro.molnr);
Nmols = numel(molnr);
Natomspermol = numel(find(gro.molnr==molnr(1)));
Ntris = size(triangles,1);

tri.verts = [gro.x gro.y gro.z];
tri.tri = repmat(triangles,[Nmols 1]) + ...
repmat(reshape(repmat((0:(Nmols-1))*Natomspermol,[Ntris 1]),[Ntris*Nmols 1]),[1 3]);
tri.n = size(tri.verts,1);

TR = triangulation(tri.tri, tri.verts);
tri.norms = vertexNormal(TR,(1:tri.n)');
