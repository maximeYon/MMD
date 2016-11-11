function ut_run
% function ut_run

rd = fullfile(fileparts(mfilename('fullpath')), '..', '..');
pd = fileparts(mfilename('fullpath'));

packages = {...
    'acq/dtd', ...
    'acq/qmas', ...
    'acq/uvec', ...
    'mdm', ...
    'mio', ...
    'models/fexi11', ...
    'models/vasco16', ...
    'models/dti_nls', ...
    'models/quick_dti'};


date_str = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
date_str = 'tmp';

log_fn = fullfile(pd, 'ut_report', sprintf('ut_%s.txt', date_str));
msf_mkdir(fileparts(log_fn));
fid = fopen(log_fn, 'w');

if (fid == -1)
    error('Could not open file for writing: %s', log_fn);
end

for c = 1:numel(packages)
    
    fprintf(fid, 'PACKAGE: %s\n', packages{c});
    fprintf(fid, '---------------------------------\n', packages{c});
    
    d = dir(fullfile(rd, packages{c}, 'ut_*.m'));
    fprintf(fid, 'Found %i ut_*.m files\n\n', numel(d));
        
    
    for c_file = 1:numel(d)
        
        m_fn = d(c_file).name(1:end-2);
        
        try
            n_ut = eval(m_fn);
        catch me
            n_ut = 0;
            disp(me.message);
        end
        
        fprintf(fid, 'Executing %i unit tests in %s\n', n_ut, m_fn);
        
        for c_ut = 1:n_ut
            try
                txt = feval(m_fn, c_ut);
                txt = [txt ': ok'];
            catch me
                txt = me.message;
            end
            
            fprintf(fid, 'Test %2.0f reports %s\n', c_ut, txt);
            
        end
        
    end
    
    fprintf(fid, '\n');
    
end

fclose(fid);