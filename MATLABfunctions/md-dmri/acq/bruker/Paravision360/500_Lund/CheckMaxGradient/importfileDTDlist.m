function LowLevelDTDr1r2list = importfileDTDlist(filename, dataLines)

%% Input handling
% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [3, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["DiffExp", "ShapeNumber", "ShapeDuration", "Gxa", "Gxb", "Gxc", "Gya", "Gyb", "Gyc", "Gza", "Gzb", "Gzc", "TE", "TR","Grel"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
LowLevelDTDr1r2list = readtable(filename, opts);

end