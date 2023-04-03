clear all

wd = cd;

datadir = '/Users/daniel/Dropbox/NMRdata/AVII500';
expnam = 'DTD'; expno = [11:121];

cd(datadir)
cd(expnam)

fidtex = fopen(['DTDfig.tex'],'w');
fprintf(fidtex,'%s\r','\documentclass{article}');
fprintf(fidtex,'%s\r','\usepackage{geometry,graphicx, epstopdf}');
fprintf(fidtex,'%s\r','\geometry{papersize={5.34in,3.3in}, scale=1}');
fprintf(fidtex,'%s\r','\usepackage[abs]{overpic}');
fprintf(fidtex,'%s\r\r','\setlength\unitlength{1in}');
fprintf(fidtex,'%s\r\r','\begin{document}');

for nexp = 1:length(expno)
    fprintf(fidtex,'%s\r','\begin{figure}');
    fprintf(fidtex,'%s\r',['  \includegraphics[width=5.34in]{' num2str(expno(nexp)) '/NII_RES/maps/dtd.pdf}']); 
    fprintf(fidtex,'%s\r\r','\end{figure}');
    fprintf(fidtex,'%s\r','\clearpage'); 
end
    
fprintf(fidtex,'%s\r','');
fprintf(fidtex,'%s\r','\end{document}');

fclose(fidtex);

