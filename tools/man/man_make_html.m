% This function makes a manual by reading the headers of all 
% matlab files in the package

cwd = fileparts(mfilename('fullpath'));
rd  = fullfile(cwd,'..','..');

packages = {...
    '.', ...
    'mdm', ...
    'mio', ...
    'mio/elastix', ...
    'models', ...
    'models/fexi11', ...
    'models/vasco16', ...
    'models/dti_nls', ...
    'models/quick_dti', ...
    'tools', ...
    'tools/tensor_maths', ...
    'tools/uvec', ...
    };

fid = fopen(fullfile(cwd, '../../index.html'),'w');
fwrite(fid, '<html><head><title>MDM manual</title></head></body><pre>');

fprintf(fid, '<h1 style="background-color:black;color:white">Packages</h1><br>');
for c = 1:numel(packages)
    if (packages{c} == '.'), continue; end
    fprintf(fid, '<a href="#p%i">%s</a><br>\n', c, packages{c});
end

for c = 1:numel(packages)
    
    fprintf(fid, '<h1 id="p%i" style="background-color:black;color:white">&nbsp;Package: %s</h1>', c, packages{c});
    
    
    % Pull in the readme file
    fid2 = fopen(fullfile(rd, packages{c}, 'readme.txt'));
    
    if (fid2 ~= -1)
        tmp = fread(fid2,inf,'char');
        warning off;
        tmp = strrep(tmp,char(10),'<br>');
        warning on;
        fprintf(fid, tmp);
        fclose(fid2);
    end
    
    fprintf(fid, '<h2 style="background-color:#dddddd">Functions</h2>');
    
    % get file headers
    d = dir(fullfile(rd, packages{c}, '*.m'));
    
    for c_file = 1:numel(d)
        
        fn2 = fullfile(rd, packages{c}, d(c_file).name);
        fid2 = fopen(fn2);
        
        % loop the m-file
        for c_line = 1:1000
            tline = fgetl(fid2);
            if ~ischar(tline), break; end
            if (numel(tline) == 0), continue; end
            
            switch (c_line)
                case 1
                    col = 'black';
                    if (~strcmpi(tline(1:min(end,8)), 'function'))
                        col = 'red';
                    end
                    
                    f_str = strtrim(tline);
                    tmp = fgetl(fid2);
                    if (~ischar(tmp)), break; end
                    
                    if (~strcmp(tmp, ['% ' f_str])), col = 'red'; end

                    fprintf(fid, '<b style="color:%s">%s</b><br>', col, tline);
                    
                    
                otherwise
                    
                    if (c_line == 2)
                        tmp = strtrim(tline);
                        if (numel(tmp) == 1) && (tmp == '%'), continue; end
                    end
                    
                    col = 'black';
                    if(tline(1) == '%')
                        fprintf(fid, '<div style="color:%s">%s</div>', col, tline);
                    else
                        break;
                    end
                    
            end
        end
        fclose(fid2);
        
        fprintf(fid, '<br>');
    end
end

fwrite(fid, '</body></html>');
