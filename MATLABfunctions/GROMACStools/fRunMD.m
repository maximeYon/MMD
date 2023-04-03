function res = fRunMD(SimPath,mdp,top,gro)

if ~isfolder(SimPath)
    res = mkdir(SimPath);
end

res = fWriteGro([SimPath '/conf.gro'],gro);
res = fWriteGro([SimPath '/restraint.gro'],gro);
res = fWriteMdp([SimPath '/grompp.mdp'],mdp);
res = fWriteTop([SimPath '/topol.top'],top);
eval(['save ' SimPath '/top.mat top'])

fnam = fullfile(SimPath,'topol.tpr');
if isfile(fnam)
    delete(fnam,[SimPath '/mdout.mdp'],[SimPath '/#*'])
end

cd(SimPath)
status = unix('source /usr/local/gromacs/bin/GMXRC; gmx grompp -quiet -maxwarn 2');

fnam = fullfile(SimPath,'traj.trr');
if isfile(fnam)
    delete(fnam,[SimPath '/md.log'],[SimPath '/ener.edr'],[SimPath '/confout.gro'],[SimPath '/#*'])
end

cd(SimPath)
status = unix('source /usr/local/gromacs/bin/GMXRC; gmx mdrun -quiet');

groout = fReadGro([SimPath '/confout.gro'],0);
res = fWriteNdx([SimPath '/index.ndx'],groout);

fnam = fullfile(SimPath,'coord.pdb');
if isfile(fnam)
    delete(fnam)
end

status = unix('source /usr/local/gromacs/bin/GMXRC; gmx traj -f confout.gro -n index.ndx -oxt coord.pdb -quiet');

%save sim.mat mdp top gro groout

res = 1;
