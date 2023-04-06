function DOR_Shape = importfile_DOR_Shape(filename, dataLines)
%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [20, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["DATA", "TYPE", "Shape", "Var4"];
opts.SelectedVariableNames = ["DATA", "TYPE", "Shape"];
opts.VariableTypes = ["double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, "Var4", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var4", "EmptyFieldRule", "auto");

% Import the data
DOR_Shape = readtable(filename, opts);

