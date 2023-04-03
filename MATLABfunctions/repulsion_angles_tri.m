clear all

Nsubdiv_v = [0];

load UDSRTriN250

wd = cd;
res = mkdir([wd '/repulsion_angles_tri']);
cd([wd '/repulsion_angles_tri'])
for n = 1:length(Nsubdiv_v);

    TR = TriRep(UDSR.tri, UDSR.x, UDSR.y, UDSR.z);
    Nsubdiv = Nsubdiv_v(n);
    TR=SubdivideSphericalMesh(TR,Nsubdiv);
    odf_s.tri = TR.Triangulation;
    odf_s.x = TR.X(:,1);
    odf_s.y = TR.X(:,2);
    odf_s.z = TR.X(:,3);
    odf_s.n = numel(odf_s.x);
    odf_s.theta = acos(odf_s.z);
    odf_s.phi = atan2(odf_s.y,odf_s.x);

    figure(1), clf
    h = trimesh(odf_s.tri,odf_s.x,odf_s.y,odf_s.z);
    set(h,'EdgeColor','b')
    axis equal
    title(['N = ' num2str(odf_s.n)])
    pause(.1)

    theta = odf_s.theta;
    phi = odf_s.phi;
    tri = odf_s.tri;
    eval(['save ' num2str(odf_s.n) ' theta phi tri'])
end
cd ..

