function paramValue=ReadPV360Param(fidPath,paramName,TwoDSeqPath)
%function paramValue=ReadPV360Param(fidPath,paramName,TwoDSeqPath)
%
%read a parameter value from the parameter files (method, acqp and reco) in
%PV 360
%
%fidPath : [string] : path to the fid file
%paramName : [string] : name of the parameter
%TwoDSeqPath : [string] : path to the 2dseq file
%
%paramValue : [string] : value found for the parameter
%
%Maxime Yon 2021

%static variables
persistent dataParsFiles;
persistent lastFidPath;
persistent lastTwoDSeqPath;

%disable missing file warnings
warning('OFF','ReadPVParam:FileNotFound');

%initialization
if(nargin==2)
    TwoDSeqPath='';
end

paramValue=[];

% if(isempty(paramName) || isempty(fidPath))
    
    %initialization
    dataParsFiles='';
    
    %path to the method and acqp files
    [dirPath]=fileparts(fidPath);
    methodPath=[dirPath,filesep,'method'];
    acqpPath=[dirPath,filesep,'acqp'];
    
    %read the method file
    fid=fopen(methodPath,'r');
    if(fid<0)
        warning('Fichier [method] introuvable !');
    else
        tmp=textscan(fid,'%s','Delimiter','\n');
        dataParsFiles=[dataParsFiles,lower(cell2mat(tmp{1}'))];
        %close the file
        fclose(fid);
    end
    
    %read the acqp file
    fid=fopen(acqpPath,'r');
    if(fid<0)
        warning('Fichier [acqp] introuvable !');
    else
        tmp=textscan(fid,'%s','Delimiter','\n');
        dataParsFiles=[dataParsFiles,lower(cell2mat(tmp{1}'))];
        %close the file
        fclose(fid);
    end
    
    if(~isempty(TwoDSeqPath))
        %path to the reco file
        [dirPath]=fileparts(TwoDSeqPath);
        recoPath=[dirPath,filesep,'reco'];
        
        %read the reco file
        fid=fopen(recoPath,'r');
        if(fid<0)
            warning('Fichier [reco] introuvable !');
        else
            tmp=textscan(fid,'%s','Delimiter','\n');
            dataParsFiles=[dataParsFiles,lower(cell2mat(tmp{1}'))];
            %close the file
            fclose(fid);
        end
    end
    
    %saving static variables
    lastFidPath=fidPath;
    lastTwoDSeqPath=TwoDSeqPath;
    
    %if we were just initializing
    if(isempty(paramName)); return; end
% end

%clean-up the parameter name
pName=['##$',lower(strtrim(paramName)),'='];

%find the position
indParam=strfind(dataParsFiles,pName);
if(isempty(indParam))
    paramValue=[];
    return;
end

%find the dimension
% pDim2=sscanf(dataParsFiles(indParam:end),[pName,'( %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d']);
EndDimVal=strfind(dataParsFiles((indParam+length(pName)):end),')');
pDimStr =dataParsFiles(indParam+length(pName)+1:indParam+length(pName)+EndDimVal(1)-2);
pDim=str2num(pDimStr);
if round(pDim)~= pDim; pDim = [];end
pDim = pDim';


%find the end of the value field
lenParamVal=strfind(dataParsFiles((indParam+length(pName)):end),'##$');
lenParamVal2=strfind(dataParsFiles((indParam+length(pName)):end),'$$ @vis');
lenParamVal=lenParamVal(1)-1;
lenParamVal2=lenParamVal2(1)-1;
lenParamVal=min(lenParamVal,lenParamVal2);
if(isempty(lenParamVal)); lenParamVal=dataParsFiles((indParam+length(pName)):end); end

%read the parameter value
if(~isempty(pDim))
    %find the beginning of the value field
    indParamVal=strfind(dataParsFiles((indParam+length(pName)):end),')');
    indParamVal=indParamVal(1);
    if(~isempty(indParamVal))
        a=indParam+length(pName)+indParamVal;
        b=a-indParamVal+lenParamVal-1;
        pValStr=dataParsFiles(a:b);
    end
else
    %find the beginning of the value field
    indParamVal=indParam+length(pName);
    a=indParamVal;
    b=a+lenParamVal-1;
    pValStr=dataParsFiles(a:b);
end

%%PV360 fix the @ values in pValStr
indAT=strfind(pValStr,'@');
while isempty(indAT)==0
    indEND = strfind(pValStr(1,indAT(1,1):end),')');
    ATparFull = pValStr(indAT(1,1):indAT(1,1)+indEND(1,1)-1);
    ATpar = pValStr(indAT(1)+1:indAT(1,1)+indEND(1)-2);
    ATpar = regexprep(ATpar,'*(',' ');
    ParNumVal=str2num(ATpar);
    ReplacePar = num2str(ones(1,ParNumVal(1,1))*ParNumVal(1,2));
    ReplacePar = [ReplacePar, ' '];
    ATparFull = regexptranslate('escape',ATparFull);
    pValStr = regexprep(pValStr,ATparFull,ReplacePar);
    indAT=strfind(pValStr,'@');
end
% End of PV360 fix the @ values in pValStr

%try converting to numeric
pValNum=str2num(pValStr);
if(isnumeric(sum(pValNum)) && ~isempty(pValNum))
    %ok: reshape
    if(~isempty(pDim))
        if length(pDim)<=2
        paramValue=squeeze(reshape(pValNum,[pDim(end:-1:1)',1]))';
        end
        if length(pDim)>2
        paramValue=permute(squeeze(reshape(pValNum,[pDim(end:-1:1)',1])),[2 1 3 4 5 6]);
        end
    else
        paramValue=pValNum;
    end
else
    %string
    paramValue=pValStr;
end

%display
if(isempty(paramValue))
    warning(['Le parametre [',paramName,'] est introuvable !']);
else
    %if numeric, concert to string
    if(isnumeric(sum(pValNum)) && ~isempty(pValNum))
        parValueForDisplay=num2str(reshape(paramValue,[1 numel(paramValue)]));
    else
        parValueForDisplay=paramValue;
    end
    
    %display
%     disp([paramName,'=',parValueForDisplay]);
end

end
