clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/opt/topspin2/data/DT/nmr';
cd(DataDir)

ExpNam = {'qMASopt'}; expno = 171:180; NDir = 8;
%ExpNam = {'isoanisocorr_qMASopt'}; expno = 27:36;

Dnocalib_array = zeros(NDir,length(expno));

cd(ExpNam{1})
for nexp = 1:length(expno)
    eval(['load ' num2str(expno(nexp)) '/Dnocalib Dnocalib'])
    Dnocalib_array(:,nexp) = Dnocalib;
end

Dmean = mean(Dnocalib_array,2);

figure(1), clf
plot(1:NDir,Dnocalib_array,'-o')
hold on
plot(1:NDir,Dmean,'k-s')

