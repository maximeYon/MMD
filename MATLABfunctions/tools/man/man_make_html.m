% This function makes a manual by reading the headers of all 
% matlab files in the package

cwd = fileparts(mfilename('fullpath'));
rd  = fullfile(cwd,'..','..');

packages = {...
    '.', ... %    'acq', ...
    'acq/bruker', ...
    'mdm', ...
    'mio', ...
    'mio/elastix', ...
    'msf', ...
    'methods', ...
    'methods/dtd',...
    'methods/dtd_pa',...
    'methods/dti_euler',...
    'methods/dti_nls',...
    'methods/erf',...
    'methods/fexi11',...
    'methods/gamma',...
    'methods/quick_dti',...
    'methods/saupe',...
    'methods/vasco16',...
    'tools', ...
    'tools/tensor_maths', ...
    'tools/uvec', ...
    };

fid = fopen(fullfile(cwd, '../../index.html'),'w');
fwrite(fid, '<html><head><title>MDM manual</title></head></body><pre>');

fprintf(fid, '<h1 style="background-color:black;color:white">Overview</h1><br>');
for c = 1:numel(packages)
    if (packages{c}(1) == '.')
        fprintf(fid, '<a href="#p%i">Introduction</a>\n', c); 
    else
        if (~exist(fullfile(rd, packages{c}), 'dir')), continue; end
        
        fprintf(fid, '<a href="#p%i">%s</a>\n', c, packages{c});
    end
end

for c = 1:numel(packages)
    
    if (~exist(fullfile(rd, packages{c}), 'dir')), continue; end
    
    fprintf(fid, '<h1 id="p%i" style="background-color:black;color:white">&nbsp;Package: %s</h1>', c, packages{c});
    
    
    % Pull in the readme file
    fid2 = fopen(fullfile(rd, packages{c}, 'readme.txt'));
    
    if (fid2 ~= -1)
        tmp = fread(fid2,inf,'char');
        warning off;
        tmp = strrep(tmp,char(10),'<br>');
        tmp = strrep(tmp, '\', '\\');
        warning on;
        fprintf(fid, char(tmp(:)'));
        fclose(fid2);
    end
    
    % get file headers
    d = dir(fullfile(rd, packages{c}, '*.m'));
    
    for c_file = 1:numel(d)
        
        if (c_file == 1)
            fprintf(fid, '<h2 style="background-color:#dddddd">Functions</h2>');
        end
        
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
                        err_msg = '1';
                    end
                    
                    % check that the name is proper
                    tmp = strtrim(tline(1:(find(tline == '(',1,'first')-1)));
                    tmp = tmp((1+max([find(tmp == '=', 1, 'last') find(tmp == ' ', 1, 'last')])):end);
                    
                    if (~strcmp(d(c_file).name, [tmp '.m']))
                        col = 'red';
                        err_msg = '2';
                    end
                    
                    f_str = strtrim(tline);
                    tmp = fgetl(fid2);
                    if (~ischar(tmp)), break; end
                    
                    if (~strcmp(tmp, ['% ' f_str])), col = 'red'; end
                    
                    if (strcmp(col, 'red'))
                        err_msg = [' (' d(c_file).name ')'];
                    else
                        err_msg = [];
                    end

                    fprintf(fid, '<b style="color:%s">%s</b>%s<br>', col, tline, err_msg);
                    
                    
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
