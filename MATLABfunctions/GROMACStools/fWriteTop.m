function res = fWriteTop(fname,top)

fid = fopen(fname,'w');

fprintf(fid,'%s\n','; Include forcefield parameters');
for count = 1:length(top.incl)
    fprintf(fid,'%s\n',['#include "' top.incl{count} '"']);
end

for nmoltyp = 1:length(top.mols)
    mol = top.mols{nmoltyp};
    if isfield(mol,'ff_custom')
        for count = 1:numel(mol.ff_custom)
            fprintf(fid,'%s\n',mol.ff_custom{count});
        end
    end
    
    fprintf(fid,'\n%s\n','[ moleculetype ]');
    fprintf(fid,'%-10s%8s\n','; Name','nrexcl');
    fprintf(fid,'%-10s%8d\n\n',top.sys.molnam{nmoltyp},3);

    fprintf(fid,'%s\n','[ atoms ]');
    fprintf(fid,'%1s%7s%10s%8s%8s%8s%8s%10s%10s\n',';','nr','type','resnr','residue','atom','cgnr','charge','mass');
    format = '%8d%10s%8d%8s%8s%8d%10.5f%10.5f\n';
    for count = 1:length(mol.atoms.nr)
        fprintf(fid,format,mol.atoms.nr(count),mol.atoms.type{count},mol.atoms.resnr(count),...
            mol.atoms.residue{count},mol.atoms.atom{count},mol.atoms.cgnr(count),...
            mol.atoms.charge(count),mol.atoms.mass(count));
    end

    fprintf(fid,'\n%s\n','[ bonds ]');
    fprintf(fid,'%1s%7s%8s%10s%12s\n',';','ai','aj','function','gromostype');
    format = '%8d%8d%10d%12s\n';
    for count = 1:length(mol.bonds.atom1)
        fprintf(fid,format,mol.bonds.atom1(count),mol.bonds.atom2(count),...
            mol.bonds.function(count),mol.bonds.gromostype{count});
    end

    fprintf(fid,'\n%s\n','[ pairs ]');
    fprintf(fid,'%1s%7s%8s%10s\n',';','ai','aj','function');
    format = '%8d%8d%10d\n';
    for count = 1:length(mol.pairs.atom1)
        fprintf(fid,format,mol.pairs.atom1(count),mol.pairs.atom2(count),...
            mol.pairs.function(count));
    end

    fprintf(fid,'\n%s\n','[ angles ]');
    fprintf(fid,'%1s%7s%8s%8s%10s%12s\n',';','ai','aj','ak','function','gromostype');
    format = '%8d%8d%8d%10d%12s\n';
    for count = 1:length(mol.angles.atom1)
        fprintf(fid,format,mol.angles.atom1(count),mol.angles.atom2(count),...
            mol.angles.atom3(count),mol.angles.function(count),mol.angles.gromostype{count});
    end

    fprintf(fid,'\n%s\n','[ dihedrals ]');
    fprintf(fid,'%1s%7s%8s%8s%8s%10s%12s\n',';','ai','aj','ak','al','function','gromostype');
    format = '%8d%8d%8d%8d%10d%12s\n';
    for count = 1:length(mol.dihedrals.atom1)
        fprintf(fid,format,mol.dihedrals.atom1(count),mol.dihedrals.atom2(count),...
            mol.dihedrals.atom3(count),mol.dihedrals.atom4(count),...
            mol.dihedrals.function(count),mol.dihedrals.gromostype{count});
    end

    fprintf(fid,'\n%s\n','[ dihedrals ]');
    fprintf(fid,'%1s%7s%8s%8s%8s%10s%12s\n',';','ai','aj','ak','al','function','gromostype');
    format = '%8d%8d%8d%8d%10d%12s\n';
    for count = 1:length(mol.impropers.atom1)
        fprintf(fid,format,mol.impropers.atom1(count),mol.impropers.atom2(count),...
            mol.impropers.atom3(count),mol.impropers.atom4(count),...
            mol.impropers.function(count),mol.impropers.gromostype{count});
    end


    fprintf(fid,'\n%s\n','#ifdef POSRES');
    fprintf(fid,'%s\n','[ position_restraints ]');
    fprintf(fid,'%1s%7s%8s%8s%8s%8s\n',';','atom','type','fx','fy','fz');
    format = '%8d%8d%8.1f%8.1f%8.1f\n';
    for count = 1:length(mol.posres.atom)
        fprintf(fid,format,mol.posres.atom(count),mol.posres.type(count),...
            mol.posres.fx(count),mol.posres.fy(count),mol.posres.fz(count));
    end
    fprintf(fid,'%s\n','#endif');
end

fprintf(fid,'\n%s\n','[ system ]');
fprintf(fid,'%s\n','; Name');
fprintf(fid,'%s\n',top.sys.name);

fprintf(fid,'\n%s\n','[ molecules ]');
fprintf(fid,'%1s%9s%8s\n',';','Compound','#mols');
for nmol = 1:length(top.sys.molnam)
    fprintf(fid,'%-10s%8d\n',top.sys.molnam{nmol},top.sys.molno{nmol});
end

fclose(fid);

res = 1;