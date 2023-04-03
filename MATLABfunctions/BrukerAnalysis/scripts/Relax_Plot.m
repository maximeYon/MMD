clear all

wd = cd;

cd(['/Users/daniel/NMRdata/AVII500/DT/'])
ExpNam = {'C10E4relax'};
expno = [91:100 51:70 27:40];
expno = [51:70];

nexp = length(expno);

cd(ExpNam{1})

for nexp = 1:length(expno)
    load([num2str(expno(nexp)) '/NMRacqus.mat']);
    load([num2str(expno(nexp)) '/RelaxDat.mat']);
    
    pl13(nexp,1) = NMRacqus.pl13;
    R(nexp,:) = FitDat.R;
    
end

expno = [71:90];

nexp = length(expno);
for nexp = 1:length(expno)
    load([num2str(expno(nexp)) '/NMRacqus.mat']);
    load([num2str(expno(nexp)) '/NutDat.mat']);
    
    pl(nexp,1) = NutDat.pl;
    B1(nexp,:) = NutDat.B1;
    
end

figure(1), clf
plref = 12;
B1ref = 38e3;

B1calc = linspace(2e3,50e3,100);
plcalc = plref + 20*log10(B1ref./B1calc);

subplot(2,2,3)
h = semilogy(pl,B1,'o',plcalc,B1calc,'-');

B1relax = B1ref*10.^(-(pl13-plref)/20);

subplot(2,2,1)
h = loglog(B1relax,R,'-o');

[Npls,Npeaks] = size(R);

subplot(2,2,2)
h = semilogy(1:Npeaks,R(1,:),'-o');

