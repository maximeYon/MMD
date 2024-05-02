function o = mdm_iter_lund(input_path, f_handle, f_args, opt)
% function o = mdm_iter_lund(input_path, f_handle, f_args, opt)

if (nargin < 1), input_path = pwd; end
if (nargin < 2), f_handle = @(a,vargin) 1; end
if (nargin < 3), f_args = {}; end
if (nargin < 4), opt = mdm_opt; end

% Assume this is the dir of individuals, filter based on some conditions
msf_log(['Searching for data in ' input_path], opt);
d = dir(input_path);

d_out = d([]);
for c = 1:numel(d)
    if any(d(c).name(1) == '._'), continue; end
    if (~d(c).isdir), continue; end
    d_out(end+1) = d(c); %#ok<AGROW>
end
d = d_out;

% Loop over the subjects (XXX: implement parfor loop)
o = cell(1, numel(d));
for c = 1:numel(d)
    o{c} = subject_loop(d(c), f_args, f_handle, input_path, opt);
end

msf_log('Finished iteration', opt);

end

% ---------------------
function o = subject_loop(d, f_args, f_handle, input_path, opt)

% Subject level
d2_dir = d.name;
d2 = dir(fullfile(input_path, d2_dir));

f_name_check = @(x) (numel(x) == 1) && (x(1) == 1);

o = {};
c_exam = 0;
for c2 = 3:numel(d2)
    
    if (~d2(c2).isdir), continue; end
    if (~f_name_check(regexp(d2(c2).name, '\d*_\d'))), continue; end
    
    msf_log(['Executing on ' d.name], opt);
    c_exam = c_exam + 1;
    
    % Exam level
    modality_list = {...
        fullfile('Diff', 'ver1'), ...
        fullfile('Diff', 'ver2'), ...
        fullfile('Diff', 'ver3'), ...
        fullfile('Diff', 'ver4'), ...
        'NII'};
    
    for c_modality = 1:numel(modality_list)
        
        d3_name = modality_list{c_modality};
        d3_dir = fullfile(d2_dir, d2(c2).name, d3_name);
        
        if (~exist(fullfile(input_path, d3_dir), 'dir')), continue, end;
        
        msf_log(['+' d3_dir], opt);
        
        s.fullfile      = fullfile(input_path, d.name, d2(c2).name, d3_name);
        s.base_path     = input_path;
        s.subject_name  = d.name;
        s.exam_name     = d2(c2).name;
        s.modality_name = d3_name;
        s.c_exam        = c_exam;
        
        try
            o{c_exam,c_modality} = f_handle(s, f_args{:});
        catch me
            msf_log(['ERROR encountered in: ' d.name ', ' d2(c2).name ', ' d3_name], opt);
            disp(getReport(me,'Extended'));
            o{c_exam,c_modality} = [];
        end
        
    end
end

end

