clear all

wd = cd;
ListDir = '/opt/topspin/exp/stan/nmr/lists/vd';

NDf = 4;
Df = logspace(log10(5e-6),log10(30e-6),NDf)
%Df = linspace(0.1e-6,1e-3,NDf);
%Df = 10e-6; NDf = length(Df);

NDm = 64;
Dm = logspace(log10(1e-6),log10(5000e-3),NDm);
%Dm = 1e-6; NDm = length(Dm);

[Df,Dm] = ndgrid(Df,Dm);
td1 = numel(Df)

Df = reshape(Df,[td1,1]);
Dm = reshape(Dm,[td1,1]);

fpath = [ListDir];

fname = 'DT_DQF_filt';
fid1 = fopen([fpath '/' fname],'w');
fname = 'DT_DQF_mix';
fid2 = fopen([fpath '/' fname],'w');

for n = 1:td1
    fprintf(fid1,'%1.8f\n',Df(n));
    fprintf(fid2,'%1.8f\n',Dm(n));
end

fclose(fid1);
fclose(fid2);
