function res = ftrr2mat(SimPath)

fnam = fullfile(SimPath,'coord.xvg');
if isfile(fnam)
    delete([SimPath '/coord.xvg'],[SimPath '/#*'])
end
cd(SimPath)
status = unix('source /usr/local/gromacs/bin/GMXRC; gmx traj -n index.ndx -ox coord.xvg -nojump -quiet');
%status = unix('g_traj -n index.ndx -ox coord.xvg -quiet');

gro = fReadGro([SimPath '/conf.gro'],0);
numperrow = (3*gro.Natom + 1);

fid = fopen([SimPath '/coord.xvg'], 'r');

format = '%f32';
datin = textscan(fid,format,'CommentStyle',{'@'},'Headerlines',13);
res = fclose(fid);

delete([SimPath '/coord.xvg'])

datin = datin{1};
Ntimes = length(datin)/numperrow;
dat = reshape(datin,numperrow,Ntimes);
clear datin

xindex = 2:3:(3*gro.Natom - 1);
yindex = 3:3:(3*gro.Natom);
zindex = 4:3:(3*gro.Natom + 1);

traj.t = dat(1,:);
traj.x = dat(xindex,:);
traj.y = dat(yindex,:);
traj.z = dat(zindex,:);
clear dat

eval(['save -v7.3 ' SimPath '/traj.mat traj'])

folder = fullfile(SimPath,'trajatom');
if ~isfolder(folder)
    res = mkdir(folder);
end

for natom = 1:gro.Natom
    trajatom.x = traj.x(natom,:);
    trajatom.y = traj.y(natom,:);
    trajatom.z = traj.z(natom,:);
    eval(['save ' SimPath '/trajatom/trajatom' num2str(natom) '.mat trajatom'])
end

res = 1;