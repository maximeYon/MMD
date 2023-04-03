function struct = fReadFFcustom(fname)

format = '%s';

fid = fopen(fname, 'r');
templine = fgetl(fid);
text = cell(1);
count = 1;
text{count,1} = templine;
while templine ~= -1
    count = count + 1;
    templine = fgetl(fid);
    text{count,1} = templine;
end

struct = text(1:(count-1));
