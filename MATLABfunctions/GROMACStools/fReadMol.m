function mol = fReadMol(MolPath)

mol.atoms = fReadAtoms([MolPath '/atoms.itp']);
mol.bonds = fReadBonds([MolPath '/bonds.itp']);
mol.pairs = fReadPairs([MolPath '/pairs.itp']);
mol.angles = fReadAngles([MolPath '/angles.itp']);
mol.dihedrals = fReadDihedrals([MolPath '/dihedrals.itp']);
mol.impropers = fReadDihedrals([MolPath '/impropers.itp']);
mol.posres = fReadPosres([MolPath '/posres.itp']);
if isfile([MolPath '/ff_custom.itp'])
    mol.ff_custom = fReadFFcustom([MolPath '/ff_custom.itp']);
end
