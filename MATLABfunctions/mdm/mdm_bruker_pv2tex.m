function res = mdm_bruker_pv2tex(data_path, tex_fn)
% function res = mdm_bruker_pv2tex(data_path, tex_fn)
%

expno = mdm_bruker_dir2expno(data_path);
%expno = expno(1:3);

subject_fn = fullfile(data_path,'subject');
fid=fopen(subject_fn,'r');
tmp=textscan(fid,'%s','Delimiter','\n');
dataParsFiles = cell2mat(tmp{1}');

expno_path = fullfile(data_path,num2str(expno(1)));

acqp = mdm_bruker_acqp2mat(expno_path);

msf_mkdir(fileparts(tex_fn));

fid = fopen(tex_fn,'w');
fprintf(fid,'%s\r','\documentclass[a4paper]{article}');
fprintf(fid,'%s\r','\usepackage{lmodern}');
fprintf(fid,'%s\r','\usepackage[T1]{fontenc}');
fprintf(fid,'%s\r','\usepackage{textcomp}');
fprintf(fid,'%s\r','\usepackage{graphicx}');
fprintf(fid,'%s\r','\usepackage{epstopdf}');
fprintf(fid,'%s\r\r','\begin{document}');

fprintf(fid,'%s\r',['\title{' acqp.user ' ' acqp.yearstr '-'...
    acqp.monthstr '-' acqp.daystr ' ' mdm_bruker_str_pv2tex(mdm_bruker_readpvparam(expno_path, 'subject_study_name')) '}']);  
ParaVal = acqp.datapath;
ParaVal(ParaVal=='_') = '-';
fprintf(fid,'%s\r\r',['\author{' ParaVal '}']);
fprintf(fid,'%s\r\r',['\maketitle ']);


for nexp = 1:length(expno)
    expno_path = fullfile(data_path,num2str(expno(nexp)));

    fprintf(fid,'%s',['\noindent \textbf{' num2str(expno(nexp)) ')} ']);

    if nexp > 1 && any(strcmp(mdm_bruker_readpvparam(expno_path, 'PULPROG'),{lower('<rFOV_DWEpiWavev1_04.ppg>'); lower('<mcw_DWEpiWavev7.ppg>')}))
        old_PVM_EpiTrajAdjkx = mdm_bruker_readpvparam(old_expno_path, 'PVM_EpiTrajAdjkx');
        PVM_EpiTrajAdjkx = mdm_bruker_readpvparam(expno_path, 'PVM_EpiTrajAdjkx');
        old_rg = mdm_bruker_readpvparam(old_expno_path, 'rg');
        rg = mdm_bruker_readpvparam(expno_path, 'rg');
        if numel(old_PVM_EpiTrajAdjkx) == numel(PVM_EpiTrajAdjkx)
            if all(old_PVM_EpiTrajAdjkx==PVM_EpiTrajAdjkx) && old_rg==rg
                fprintf(fid,'%s','\textbf{GOP}, ');
            end
        end        
    end
    
    ParaToList = {'pulprog','EchoTime','DwLoopOrder','DwNAmplitudes','DwGradShapeStrArr','DwGradShapeStrArr1','DwGradShapeStrArr2','DwNDirs','DwGradDur','DwGradDur1','DwGradDur2','DwGradTsep','rg'};
    for npara = 1:length(ParaToList)
        ParaVal = mdm_bruker_readpvparam(expno_path, ParaToList{npara});
        if exist('old_expno_path')
            oldParaVal = mdm_bruker_readpvparam(old_expno_path, ParaToList{npara});
        else
            oldParaVal = eps+.001;
        end
        
        if ischar(ParaVal) 
            if ~strcmp(ParaVal,oldParaVal) && ~isempty(ParaVal)
                fprintf(fid,'%s',[mdm_bruker_str_pv2tex(ParaToList{npara}) ' = ' mdm_bruker_str_pv2tex(ParaVal) ', ']);
            end
        else
            if ~isempty(ParaVal)
                if isempty(oldParaVal), oldParaVal = eps+.001; end
                if ParaVal~=oldParaVal
                    fprintf(fid,'%s',[mdm_bruker_str_pv2tex(ParaToList{npara}) ' = ' num2str(ParaVal) ', ']);
                end
            end
        end
        
    end
        
    acqp = mdm_bruker_acqp2mat(expno_path);
    fprintf(fid,'%s',[acqp.hourstr ':' acqp.minutestr]);

    fprintf(fid,'.\r\r');

    old_expno_path = expno_path;
end

fprintf(fid,'%s\r','');
fprintf(fid,'%s\r','\end{document}');

fclose(fid);

res = 1;
