function res = fPlotCylindersTriMesh(gro,cylpars,figno)

if nargin < 3
    figno = 1;
end

% gro = fPutMoleculesInBoxGro(gro);
% gro = fAlignGyrationTensorGro(gro);
% gro = fCOMToOriginGro(gro);        
% %gro.y = gro.y - min(gro.y);

tri = fCalcTriSurfParsCylindersTriMesh(gro,cylpars);

TR = triangulation(tri.tri, tri.verts);

figure(figno), clf
trisurf(TR,'FaceColor',.9*[1 1 1],'EdgeColor','k');
axis equal
view(30,30)
xlabel('x'), ylabel('y'), zlabel('z')

res = 1;