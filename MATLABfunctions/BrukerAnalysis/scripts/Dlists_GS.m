clear all

wd = cd;
ListDir = '/opt/topspin/exp/stan/nmr/lists/vd';

NDf = 8;
Df = logspace(-7,log10(1e-3),NDf);
%Df = linspace(0.1e-6,1e-3,NDf);

NDm = 32;
Dm = logspace(-4,log10(5),NDm);

[Df,Dm] = ndgrid(Df,Dm);
td1 = numel(Df)

Df = reshape(Df,[td1,1]);
Dm = reshape(Dm,[td1,1]);

fpath = [ListDir];

fname = 'DT_GS_filt';
fid1 = fopen([fpath '/' fname],'w');
fname = 'DT_GS_mix';
fid2 = fopen([fpath '/' fname],'w');

for n = 1:td1
    fprintf(fid1,'%1.8f\n',Df(n));
    fprintf(fid2,'%1.8f\n',Dm(n));
end

fclose(fid1);
fclose(fid2);
