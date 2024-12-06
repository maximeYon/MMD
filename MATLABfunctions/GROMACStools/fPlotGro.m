function res = fPlotGro(gro,mol,figno)

Nmol = max(gro.molnr);
index = (1:length(mol.atoms.nr))';
molnr = 1;

figure(figno), clf
for nmol = 1:Nmol
    Plot.Atoms.x = gro.x(index);
    Plot.Atoms.y = gro.y(index);
    Plot.Atoms.z = gro.z(index);

    Plot.Bonds.x = [Plot.Atoms.x(mol.bonds.atom1)'; Plot.Atoms.x(mol.bonds.atom2)'];
    Plot.Bonds.y = [Plot.Atoms.y(mol.bonds.atom1)'; Plot.Atoms.y(mol.bonds.atom2)'];
    Plot.Bonds.z = [Plot.Atoms.z(mol.bonds.atom1)'; Plot.Atoms.z(mol.bonds.atom2)'];

    plot3(Plot.Bonds.x,Plot.Bonds.y,Plot.Bonds.z,'k-')
    hold on
    if gro.Natom < 20
        plot3(Plot.Atoms.x,Plot.Atoms.y,Plot.Atoms.z,'ko')
    end
    axis equal
    
    index = index + length(mol.atoms.nr);
    molnr = molnr+1;
end

axis([0 gro.boxxx 0 gro.boxyy 0 gro.boxzz])

res = 1;