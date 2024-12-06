clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/Documents/Spaces/Presentations';

ExpNam = {'AOToct_Eq3'}; expno = 15:10:185;
ExpNam = {'AOToct_Eq4'}; expno = [15:10:945];
ExpNam = {'AOToct_Eq5'}; expno = [5:10:485];
ExpNam = {'AOToct_Eq6'}; expno = [15:10:245];
%ExpNam = {'AOToct_Eq7'}; expno = [15:10:245]; %expno = 25;

cd(DataDir)
cd(ExpNam{1})

fidtex = fopen(['ODFFig.tex'],'w');
fprintf(fidtex,'%s\r','\documentclass{article}');
fprintf(fidtex,'%s\r','\usepackage{geometry,graphicx, epstopdf}');
fprintf(fidtex,'%s\r','\geometry{papersize={20in,20in}, scale=1}');
fprintf(fidtex,'%s\r','\usepackage[abs]{overpic}');
fprintf(fidtex,'%s\r\r','\setlength\unitlength{1in}');
fprintf(fidtex,'%s\r\r','\begin{document}');

for nexp = 1:length(expno)
    fprintf(fidtex,'%s\r','\begin{figure}');
    fprintf(fidtex,'%s\r',['  \includegraphics[width=20in]{' num2str(expno(nexp)) '/ODFFig.pdf}']); 
    fprintf(fidtex,'%s\r\r','\end{figure}');
    fprintf(fidtex,'%s\r','\clearpage'); 
end
    
fprintf(fidtex,'%s\r','');
fprintf(fidtex,'%s\r','\end{document}');

fclose(fidtex);

