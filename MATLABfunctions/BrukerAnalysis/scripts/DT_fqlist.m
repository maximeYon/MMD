clear all

wd = cd;

DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';

fqmin = -70e3; fqmax = 70e3; Nfq = 15;
%fqmin = 80e3; fqmax = -70e3; Nfq = 16;
%fqmin = -30e3; fqmax = 30e3; Nfq = 7;
%fqmin = -120e3; fqmax = 120e3; Nfq = 25;

fq = linspace(fqmin,fqmax,Nfq);


ListDir = '/opt/topspin2/exp/stan/nmr/lists/f1';

fpath = [ListDir];

fname = 'DT_fq';
fid = fopen([fpath '/' fname],'w');

fprintf(fid,'%s\n','sfo hz');
for n = 1:Nfq
    fprintf(fid,'%1.2f\n',fq(n));
end

fclose(fid);


cd(wd)