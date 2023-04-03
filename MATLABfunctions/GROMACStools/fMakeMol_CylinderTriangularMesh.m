function [gro,mol] = fMakeMol_CylinderTriangularMesh(InVars,MolPath)

MolNam = InVars.MolNam;
Nphi = InVars.Nphi;
Nz = InVars.Nz;
bondlength = InVars.bondlength;

atommass = InVars.atommass;
bond_springconstant = InVars.bond_springconstant;
angle_springconstant = InVars.angle_springconstant;
rm = InVars.LJ_well_distance; % distance at potential well
epsilon = InVars.LJ_well_depth; % depth of potential well


c6 = 2*epsilon*rm^6;
c12 = epsilon*rm^12;

cylinderradius = bondlength/(2*sin(pi/Nphi));
cylinderlength = sqrt(3)/2*Nz*bondlength;
aspectratio = cylinderlength/2/cylinderradius;

bondangle_radial = (Nphi-2)*pi/Nphi;
ri = bondlength/2*[cos(pi/Nphi)/sin(pi/Nphi) -1 -sqrt(3)];
rj = bondlength/2*[           1/sin(pi/Nphi)  0        0];
rk = bondlength/2*[cos(pi/Nphi)/sin(pi/Nphi)  1  sqrt(3)];
rij = rj-ri;
rjk = rk-rj;
bondangle_diagonal = pi - acos(sum(rij.*rjk)/(sqrt(sum(rij.*rij))*sqrt(sum(rjk.*rjk))));

phi = 2*pi*linspace(0,1,Nphi+1)';
dphi = phi(2) - phi(1);
phi = phi(1:Nphi);

z = cylinderlength*linspace(0,1,Nz+1)';
z = z(1:Nz);

[phi,z] = ndgrid(phi,z);
[phi_count,z_count] = ndgrid(1:Nphi,1:Nz);
phi = phi + (z_count-1)*dphi/2;

x = cylinderradius*cos(phi);
y = cylinderradius*sin(phi);

x = x(:); y = y(:); z = z(:);

Natoms = numel(x);

bonds = [];
for nz = 1:Nz
    bonds_temp = [(1:Nphi)' circshift((1:Nphi)',-1,1)] + (nz-1)*Nphi;
    bonds = cat(1,bonds,bonds_temp);
end
for nphi = 1:Nphi
    atoms_temp = (1+(0:Nz-1)*Nphi)' + (nphi-1);
    bonds_temp = [atoms_temp(1:(end-1)) atoms_temp(2:end)];
    bonds = cat(1,bonds,bonds_temp);
end
for nphi = 1:Nphi
    atoms_temp = 1 + (((1:Nz)-1)*(Nphi-1))' + (nphi-1);    
    for ncycle = 1:ceil(Nz/Nphi)
        ind = (nphi+1 + (ncycle-1)*Nphi):Nz;
        atoms_temp(ind) = atoms_temp(ind) + Nphi;
    end    
    bonds_temp = [atoms_temp(1:(end-1)) atoms_temp(2:end)];
    bonds = cat(1,bonds,bonds_temp);
end
Nbonds = size(bonds,1);


angles = [];
for nz = 1:Nz
    angles_temp = [(1:Nphi)' circshift((1:Nphi)',-1,1) circshift((1:Nphi)',-2,1)] + (nz-1)*Nphi;
    angles = cat(1,angles,angles_temp);
end
Nangles_radial = size(angles,1);
for nphi = 1:Nphi
    atoms_temp = (1+(0:Nz-1)*Nphi)' + (nphi-1);
    angles_temp = [atoms_temp(1:(end-2)) atoms_temp(2:(end-1)) atoms_temp(3:end)];
    angles = cat(1,angles,angles_temp);
end
for nphi = 1:Nphi
    atoms_temp = 1 + (((1:Nz)-1)*(Nphi-1))' + (nphi-1);    
    for ncycle = 1:ceil(Nz/Nphi)
        ind = (nphi+1 + (ncycle-1)*Nphi):Nz;
        atoms_temp(ind) = atoms_temp(ind) + Nphi;
    end    
    angles_temp = [atoms_temp(1:(end-2)) atoms_temp(2:(end-1)) atoms_temp(3:end)];
    angles = cat(1,angles,angles_temp);
end
Nangles = size(angles,1);
Nangles_diagonal = Nangles - Nangles_radial;

triangles = [];
for nz = 1:(Nz-1)
    triangles_temp = [(1:Nphi)' circshift((1:Nphi)',-1,1) (1:Nphi)'+Nphi] + (nz-1)*Nphi;
    triangles = cat(1,triangles,triangles_temp);
end
for nz = 1:(Nz-1)
    triangles_temp = [(1:Nphi)' (1:Nphi)'+Nphi circshift((1:Nphi)'+Nphi,1,1) ] + (nz-1)*Nphi;
    triangles = cat(1,triangles,triangles_temp);
end
Ntriangles = size(triangles,1);


% figure(1), clf
% plot3(x,y,z,'.')
% axis equal
% hold on
% plot3(x(1),y(1),z(1),'x')
% for nbond = 1:Nbonds
%     plot3(x(bonds(nbond,:)),y(bonds(nbond,:)),z(bonds(nbond,:)),'k-')
% end
% 
% for nangle = 1:Nangles
%     plot3(x(angles(nangle,:)),y(angles(nangle,:)),z(angles(nangle,:)),'r-')
% end
% 
% for ntriangle = 1:Ntriangles
%     plot3(x(triangles(ntriangle,:)),y(triangles(ntriangle,:)),z(triangles(ntriangle,:)),'g-')
% end

if ~isdir(MolPath)
    mkdir(MolPath)
end

mol.atoms.nr = zeros(Natoms,1);
mol.atoms.type = cell(Natoms,1); 
mol.atoms.resnr = zeros(Natoms,1);
mol.atoms.residue = cell(Natoms,1); 
mol.atoms.atom = cell(Natoms,1); 
mol.atoms.cgnr = zeros(Natoms,1);
mol.atoms.charge = zeros(Natoms,1);
mol.atoms.mass = zeros(Natoms,1);
for natom = 1:Natoms
    mol.atoms.nr(natom,1) = natom;
    mol.atoms.type{natom,1} = 'AX';
    mol.atoms.resnr(natom,1) = 1;
    mol.atoms.residue{natom,1} = MolNam;
    mol.atoms.atom{natom,1} = ['C' num2str(natom)];
    mol.atoms.cgnr(natom,1) = ceil(natom/4);
    mol.atoms.charge(natom,1) = 0;
    mol.atoms.mass(natom,1) = atommass;
end

fid = fopen(fullfile(MolPath,'atoms.itp'),'w');
fprintf(fid,'%1s%7s%10s%8s%8s%8s%8s%10s%10s\n',';','nr','type','resnr','residue','atom','cgnr','charge','mass');
format = '%8d%10s%8d%8s%8s%8d%10.5f%10.5f\n';
for count = 1:length(mol.atoms.nr)
    fprintf(fid,format,mol.atoms.nr(count),mol.atoms.type{count},mol.atoms.resnr(count),...
        mol.atoms.residue{count},mol.atoms.atom{count},mol.atoms.cgnr(count),...
        mol.atoms.charge(count),mol.atoms.mass(count));
end
fclose(fid);

mol.bonds.atom1 = bonds(:,1);
mol.bonds.atom2 = bonds(:,2);
mol.bonds.function = 2*ones(Nbonds,1);
mol.bonds.gromostype = cell(Nbonds,1); 
for nbond = 1:Nbonds
    mol.bonds.gromostype{nbond,1} = 'gb_rad';
end

fid = fopen(fullfile(MolPath,'bonds.itp'),'w');
fprintf(fid,'%1s%7s%8s%10s%12s\n',';','ai','aj','function','gromostype');
format = '%8d%8d%10d%12s\n';
for count = 1:length(mol.bonds.atom1)
    fprintf(fid,format,mol.bonds.atom1(count),mol.bonds.atom2(count),...
        mol.bonds.function(count),mol.bonds.gromostype{count});
end
fclose(fid);

mol.angles.atom1 = angles(:,1);
mol.angles.atom2 = angles(:,2);
mol.angles.atom3 = angles(:,3);
mol.angles.function = 2*ones(Nangles,1);
mol.angles.gromostype = cell(Nangles,1); 
for nangle = 1:Nangles_radial
    mol.angles.gromostype{nangle,1} = 'ga_rad';
end
for nangle = (Nangles_radial+1):Nangles
    mol.angles.gromostype{nangle,1} = 'ga_diag';
end

fid = fopen(fullfile(MolPath,'angles.itp'),'w');
fprintf(fid,'%1s%7s%8s%8s%10s%12s\n',';','ai','aj','ak','function','gromostype');
format = '%8d%8d%8d%10d%12s\n';
for count = 1:length(mol.angles.atom1)
    fprintf(fid,format,mol.angles.atom1(count),mol.angles.atom2(count),...
        mol.angles.atom3(count),mol.angles.function(count),mol.angles.gromostype{count});
end
fclose(fid);

mol.pairs.atom1 = [];
mol.dihedrals.atom1 = [];
mol.impropers.atom1 = [];

fid = fopen(fullfile(MolPath,'pairs.itp'),'w');
fprintf(fid,'%1s%7s%8s%10s\n',';','ai','aj','function');
format = '%8d%8d%10d\n';
for count = 1:length(mol.pairs.atom1)
    fprintf(fid,format,mol.pairs.atom1(count),mol.pairs.atom2(count),...
        mol.pairs.function(count));
end
fclose(fid);

fid = fopen(fullfile(MolPath,'dihedrals.itp'),'w');
fprintf(fid,'%1s%7s%8s%8s%8s%10s%12s\n',';','ai','aj','ak','al','function','gromostype');
format = '%8d%8d%8d%8d%10d%12s\n';
for count = 1:length(mol.dihedrals.atom1)
    fprintf(fid,format,mol.dihedrals.atom1(count),mol.dihedrals.atom2(count),...
        mol.dihedrals.atom3(count),mol.dihedrals.atom4(count),...
        mol.dihedrals.function(count),mol.dihedrals.gromostype{count});
end
fclose(fid);

fid = fopen(fullfile(MolPath,'impropers.itp'),'w');
fprintf(fid,'%1s%7s%8s%8s%8s%10s%12s\n',';','ai','aj','ak','al','function','gromostype');
format = '%8d%8d%8d%8d%10d%12s\n';
for count = 1:length(mol.impropers.atom1)
    fprintf(fid,format,mol.impropers.atom1(count),mol.impropers.atom2(count),...
        mol.impropers.atom3(count),mol.impropers.atom4(count),...
        mol.impropers.function(count),mol.impropers.gromostype{count});
end
fclose(fid);

mol.posres.atom = 1;
mol.posres.type = 1;
mol.posres.fx = 0;
mol.posres.fy = 0;
mol.posres.fz = 1000;

fid = fopen(fullfile(MolPath,'posres.itp'),'w');
fprintf(fid,'%1s%7s%8s%8s%8s%8s\n',';','atom','type','fx','fy','fz');
format = '%8d%8d%8.1f%8.1f%8.1f\n';
for count = 1:length(mol.posres.atom)
    fprintf(fid,format,mol.posres.atom(count),mol.posres.type(count),...
        mol.posres.fx(count),mol.posres.fy(count),mol.posres.fz(count));
end
fclose(fid);

ff_custom = {['#define gb_rad ' num2str(bondlength,4) ' ' num2str(bond_springconstant)];
    ['#define ga_rad ' num2str(bondangle_radial/pi*180,4) ' ' num2str(angle_springconstant)];
    ['#define ga_diag ' num2str(bondangle_diagonal/pi*180,4) ' ' num2str(angle_springconstant)];
    ['[ atomtypes ]'];
    [';name  at.num   mass      charge  ptype       c6           c12'];
    ['  AX    6 	0.000      0.000     A ' num2str(c6,6) ' ' num2str(c12,6)]};

fid = fopen(fullfile(MolPath,'ff_custom.itp'),'w');
fprintf(fid,'%s\n','; Custom forcefield parameters');
for count = 1:size(ff_custom,1)
    fprintf(fid,'%s\n',ff_custom{count});
end
fclose(fid);

gro.header = [MolNam ' starting configuration'];
gro.Natom = Natoms;
gro.molnr = zeros(Natoms,1);
gro.molecule = cell(Natoms,1);
gro.atom = cell(Natoms,1);
gro.atomnr = zeros(Natoms,1);

gro.boxxx = max([4*cylinderradius 3]);
gro.boxyy = max([4*cylinderradius 3]);
gro.boxzz = max([2*cylinderlength 3]);

gro.x = zeros(Natoms,1);
gro.y = zeros(Natoms,1);
gro.z = zeros(Natoms,1);

for natom = 1:Natoms
    gro.molnr(natom,1) = 1;
    gro.molecule{natom,1} = MolNam;
    gro.atom{natom,1} = ['C' num2str(natom)];
    gro.atomnr(natom,1) = natom;
    gro.x(natom,1) = x(natom) + .5*gro.boxxx;
    gro.y(natom,1) = y(natom) + .5*gro.boxyy;
    gro.z(natom,1) = z(natom) + 0.25*gro.boxzz;
end

fWriteGro(fullfile(MolPath,'conf.gro'),gro);

return
membrane.verts = [x y z];
membrane.tri = triangles;
membrane.n = Natoms;

TR = triangulation(membrane.tri, membrane.verts);
membrane.norms = vertexNormal(TR,(1:membrane.n)');

trisurf(TR,'FaceColor',[0.8 0.8 1.0]);
axis equal
hold on
quiver3(membrane.verts(:,1),membrane.verts(:,2),membrane.verts(:,3), ...
     membrane.norms(:,1),membrane.norms(:,2),membrane.norms(:,3),0.5,'Color','b');
 
bondlength
cylinderradius
cylinderlength
numel(z)


clear tri
tri.verts = membrane.verts;
tri.tri = membrane.tri;

wd = pwd;

povinc_path = fullfile(wd,'tri');
pov_tri2povinc(tri,povinc_path);

