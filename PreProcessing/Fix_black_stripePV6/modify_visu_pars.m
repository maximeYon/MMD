function []=modify_visu_pars(Path,paramName,new_values)

filename=[Path,filesep,'visu_pars'];

%read the method file
fileID=fopen(filename,'r');
[~,~,machinefmt,~] = fopen(fileID);
if(fileID<0)
    warning('Fichier visu_pars introuvable !');
else
    formatSpec = '%s%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
    %% Remove white space around all cell columns.
    dataArray{1} = dataArray{1};
    fclose(fileID);
end
clearvars fileID;

%clean-up the parameter name
pName=['##$',strtrim(paramName),'='];

%find the position
pos_ind = zeros(1,size(dataArray{1,1},1));
for ind = 1:size(dataArray{1,1},1)
    pos_ind(1,ind) = numel(strfind(dataArray{1, 1}{ind, 1}, pName));
end
pos_param = find(pos_ind);

%find the position of the next parameter
pos_ind = zeros(1,size(dataArray{1,1},1));
for ind = pos_param+1:size(dataArray{1,1},1)
    pos_ind(1,ind) = numel(strfind(dataArray{1, 1}{ind, 1}, '##$'));
end
pos_next_param = find(pos_ind);
pos_next_param = pos_next_param(1,1);
size_param = pos_next_param-pos_param-1;

%% Write in array
if numel(strfind(dataArray{1, 1}{pos_param, 1}, '=('))
    if size_param ~=1
        new_values = reshape(new_values,size(new_values,2)/size_param,size_param);
        ind = 0;
        my_form = '%.18f ';
        my_form2 = '%.18f ';
        for indform =1:size(new_values,1)-1
            my_form = [my_form my_form2];
        end
        for ind_p = pos_param +1:pos_param+size_param
            ind = ind+1;
            dataArray{1, 1}{ind_p, 1} = [num2str(new_values(:,ind)',my_form) ' '];
        end
    else
    dataArray{1, 1}{pos_param +1, 1} = num2str(new_values,'%d ');
    end
else
    dataArray{1, 1}{pos_param, 1} = [pName num2str(new_values,'%d ')];
end

%% Write new file
% warning off MATLAB:iofun:UnsupportedEncoding;
% fileID = fopen(filename, 'w', 'b', 'UTF16-LE');
fileID=fopen(filename,'w',machinefmt);
formatSpec = '%s\n';
for ind = 1:size(dataArray{1,1},1)
    fprintf(fileID,formatSpec,dataArray{1, 1}{ind, 1});
end
fclose(fileID);
end
