function res = fWriteMdp(fname,mdp)

fid = fopen(fname,'w');

names = fieldnames(mdp);

for count = 1:length(names)
    eval(['paranam = mdp.' names{count} '{1};'])
    eval(['paraval = num2str(mdp.' names{count} '{2});'])
    fprintf(fid,'%-25s%s%s\n',paranam,' = ',paraval);
end

fclose(fid);

res = 1;