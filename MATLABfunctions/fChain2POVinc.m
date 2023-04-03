function res = fChain2POVinc(fnam,Chain)

fid = fopen(fnam,'w');

Nnodes = length(Chain.x);
tempstr = ['sphere{<' num2str(Chain.x(1)) ','...
    num2str(Chain.z(1)) ',' ...
    num2str(Chain.y(1)) '>' ',rad}']; 
fprintf(fid,'%s\n',tempstr);
for nnode = 2:Nnodes
    tempstr = ['cylinder{<' num2str(Chain.x(nnode-1)) ','...
        num2str(Chain.z(nnode-1)) ',' ...
        num2str(Chain.y(nnode-1)) '>,<'...
        num2str(Chain.x(nnode)) ','...
        num2str(Chain.z(nnode)) ',' ...
        num2str(Chain.y(nnode)) '>,rad}']; 
    fprintf(fid,'%s\n',tempstr);
    tempstr = ['sphere{<' num2str(Chain.x(nnode)) ','...
        num2str(Chain.z(nnode)) ',' ...
        num2str(Chain.y(nnode)) '>' ',rad}']; 
    fprintf(fid,'%s\n',tempstr);
end

fclose(fid);

res = 1;