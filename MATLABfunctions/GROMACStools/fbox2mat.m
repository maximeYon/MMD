function res = fbox2mat(SimPath)

delete([SimPath '/box.xvg'],[SimPath '/#*'])
cd(SimPath)

numperrow = 7;

status = unix('source /usr/local/gromacs/bin/GMXRC; gmx traj -n index.ndx -ob box.xvg -quiet');

fid = fopen([SimPath '/box.xvg'], 'r');
format = '%f32';

datin = textscan(fid,format,'CommentStyle',{'@'},'Headerlines',13);
datin = datin{1};
Ntimes = length(datin)/numperrow;
dat = reshape(datin,numperrow,Ntimes);

box.t = dat(1,:);
box.xx = dat(2,:);
box.yy = dat(3,:);
box.zz = dat(4,:);
box.yx = dat(5,:);
box.zx = dat(6,:);
box.zy = dat(7,:);

eval(['save ' SimPath '/box.mat box'])

res = fclose(fid);

delete([SimPath '/box.xvg'])

res = 1;