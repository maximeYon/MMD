function mdp = fReadMdp(fname)

fid = fopen(fname, 'r');

structnam = 'mdp';

templine = fgetl(fid);
while ischar(templine) 
    if ~isempty(templine)
        if ~strcmp(';',templine(1))
            
            indxeq = find(templine == '=');
            
            paranam = deblank(templine(1:indxeq(1)-1));
            paranammdp = paranam;
            indxdash = find(paranam == '-');
            if ~isempty(indxdash)
                paranam(indxdash) = '_';
            end
            
            tempval = templine(indxeq+1:length(templine));
            tempval = deblank(strjust(tempval, 'left'));
            indxsemicol = find(tempval == ';');
            if ~isempty(indxsemicol)
                tempval = tempval(1:(min(indxsemicol)-1));
            end
            tempval = deblank(strjust(tempval, 'left'));
            tf = isstrprop(tempval, 'alpha');
            indxalpha = find(tf==1);
            tempval1 = tempval;
            if all(tf == 0);
                tempval = str2num(tempval);
            elseif all(tempval(indxalpha)=='e')
                tempval = str2num(tempval);
            end
            if strcmp(tempval1,'9-3')
                tempval = tempval1;
            end

            
            eval([structnam '.' paranam ' = {paranammdp, tempval};']);
        end
    end
    
    templine = fgetl(fid);
end
 
fclose(fid);