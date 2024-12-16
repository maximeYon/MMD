function status = mdm_fit(varargin)
% function status = mdm_fit(varargin)
%
% run function without argument to see help text

% initiate
status = 0; % 0 - error, 1 - ok
nargin = numel(varargin);

if (nargin == 0), print_msg('No arguments supplied'); return; end

input.present = 1;

mandatory_arguments = {'data', 'out', 'method'};
optional_arguments  = {'mask', 'gwf', 'btens', 'xps', 'verbose', ...
    'i_range', 'j_range', 'k_range'};

% parse input arguments
c = 1;
nargin = numel(varargin);
while (c <= nargin)
    
    if (strcmp(varargin{c}, '--help')) % show help msg, terminate
        print_msg('Help requested'); return;
    end
    
    if (~strcmp(varargin{c}(1:min(end,2)), '--')) % specify input
        print_msg('Incorrectly specified arguments'); return;
    end
        
    input_flag = varargin{c}(3:end); c = c + 1;
    
    if (c > nargin)
        print_msg(['Insufficient number of arguments, ' input_flag]); return;
    end        
    
    data_value = varargin{c}; c = c + 1;    
    
    
    if (~strcmpi(input_flag, cat(2, mandatory_arguments, optional_arguments)))
        print_msg(['Unknown argument: ' input_flag]); return;
    end
    
    switch (input_flag)
        case {'i_range', 'j_range', 'k_range'}
            if (ischar(data_value))
                data_value = str2num(data_value);
            end
    end
            
    
    
    input.(input_flag) = data_value;
    
end

% check that all mandatory arguments are present
for c = 1:numel(mandatory_arguments)
    if (~isfield(input, mandatory_arguments{c}))
        print_msg(['Missing argument: ' mandatory_arguments{c}]); return;
    end
end

% check that we get either a b-tensor file, a gwf-set, or an xps
metadata_input = [isfield(input, 'btens') ...
    isfield(input, 'gwf') isfield(input, 'xps')];

if (sum(metadata_input) ~= 1)
    print_msg('Missing input (specify either btens, gwf or xps)'); return;
end    


% ---------------------------------------------------------------------
% build xps
% ---------------------------------------------------------------------
clear xps;
if (isfield(input, 'xps'))
    
    try
        xps = mdm_xps_load(input.xps);
    catch me
        print_msg(['Could not load xps file (' me.message ')']); return;
    end
    
elseif (isfield(input, 'gwf'))
    
        print_msg('Not yet implemented'); return;
        
elseif (isfield(input, 'btens'))
    
    try
        xps = mdm_xps_from_btens(input.btens);
    catch me    
        print_msg(['Could not load btens file (' me.message ')']); return;
    end        
end

% ------------------------------------------------------------------------
% Convert input structure to s-structure and xps
% ------------------------------------------------------------------------

% Ensure that input file exists
s.nii_fn = input.data;
s.xps = xps;
clear xps;
if (~exist(s.nii_fn, 'file'))
    print_msg(['Data file not found (' s.nii_fn ')']); return;
end

% Check filetype
[~,~,ext] = msf_fileparts(s.nii_fn);
if (~any(strcmpi(ext, {'.nii', '.nii.gz'})))
    print_msg(['Unknown format: ' ext ', only .nii/.nii.gz allowed']); return;
end

% Connect to mask
if (isfield(input, 'mask'))
    s.mask_fn = input.mask;
    
    if (~exist(s.mask_fn, 'file'))
        print_msg(['Mask not found (' s.mask_fn ')']); return;
    end
    
    % Check filetype
    [~,~,ext] = msf_fileparts(s.mask_fn);
    if (~any(strcmpi(ext, {'.nii', '.nii.gz'})))
        print_msg(['Unknown format o mask: ' ext ', only .nii/.nii.gz allowed']);
        return;
    end
end


% check if method is present and that the constructed xps is sufficient
m = input.method;
m(~ismember(m, ['a':'z' 'A':'Z' '0':'9' '_'])) = ''; % protect

if (~exist([m '_check_xps'], 'file'))
    print_msg(['Did not find method ' input.method]);
    return;
end    

try
    feval([m '_check_xps'], s.xps);
catch me
    print_msg(['Check input data (' me.message ')']); return;
end



% work with the output structure
[out_path, out_prefix] = fileparts(input.out);
if (~isempty(out_prefix) && (out_prefix(end) ~= '_')) 
    out_prefix = [out_prefix '_'];
end

% ------------------------------------------------------------------------
% Start analysis
% ------------------------------------------------------------------------
try
    opt = feval([m '_opt'], mdm_opt());
    
    % control fit region (for debugging/quick testing)
    if (isfield(input,'i_range')), opt.i_range = input.i_range; end
    if (isfield(input,'j_range')), opt.j_range = input.j_range; end
    if (isfield(input,'k_range')), opt.k_range = input.k_range; end
    
    opt.do_overwrite = 1;
    
    paths = mdm_paths(out_path, out_prefix, ['_' m]);
    
    % If needed, run powder averaging
    if (isfield(opt.(m), 'pa_method')) && (opt.(m).pa_method)
        s = mdm_s_powder_average(s, paths.nii_path, opt);
    end
    
    % Run the analysis
    mdm_data2fit(eval(['@' m '_4d_data2fit']), s, paths.mfs_fn, opt);
    mdm_fit2param(eval(['@' m '_4d_fit2param']), paths.mfs_fn, paths.dps_fn, opt);
    
    % Save nifti parameter maps
    if (~isfield(opt.(m),'fig_prefix')), 
        opt.(m).fig_prefix = [m '_'];
    end
    
    opt.(m).fig_prefix = [out_prefix opt.(m).fig_prefix];
    mdm_param2nii(paths.dps_fn, paths.nii_path, opt.(m), opt);
    
catch me
    print_err_message(me); return; 
end

% Everything finished OK
status = 1;

end


% -------------------------------------------------------------------------
% Quit the program with an error message and show the help text
% -------------------------------------------------------------------------
function print_msg(msg)

if (nargin < 1), msg = ''; end

msg = strrep(msg, '\', '\\');

msg_correct_usage = [...
    'Command incorrectly used: ' msg '\n' ...
    '\n' ...
    'Correct usage: mdm_fit FLAGS\n' ...
    '\n' ...
    '  where FLAGS is\n' ...
    '\n' ...
    '     --data filename (to nifti volume)\n' ...
    '\n' ...
    '     --method name (one of methods below)\n' ...
    '          dtd_gamma           DIVIDE approach\n' ...
    '          dtd_codivide        CODIVIDE approach\n' ...
    '          dtd_ndi             Simplified NODDI\n' ...
    '          fexi11              Filter exchange imaging\n' ...
    '          ...                 See methods/* for more\n' ...
    '\n' ...
    '     --out output (path and prefix as path/prefix)\n' ...
    '\n' ....
    '     --btens/xps/gwf file (use one of the three)\n' ...
    '          btens               b-tensor file in voigt format (unit: um2/ms)\n' ...
    '          xps                 xps file as .mat\n' ...
    '          gwf                 gradient waveforms (not implemented)\n' ...
    '\n' ...
    'Optional: \n' ...
    '\n' ...
    '     -mask filename           Mask for the computations\n' ...
    '\n' ...
    'Other: \n' ...
    '\n' ...
    '     -help                    Display help\n' ...
    ];

disp(sprintf(msg_correct_usage));

end

% -------------------------------------------------------------------------
% Print the error message for easier debugging
% -------------------------------------------------------------------------
function print_err_message(me)

disp(['Error: ' me.message ' (' me.stack(1).file ')']);

disp('Stack:');
for c = 1:numel(me.stack)
    str = [me.stack(c).name ' at ' num2str(me.stack(c).line)];
    disp(str);
end

end
