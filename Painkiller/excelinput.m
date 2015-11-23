
% Convert table to cell array of structs
% Columns that begin with nonlin are concatenated to a single string
specs = readtable('testinput.xlsx')
specs2 = readtable('testinput.xlsx', 'Sheet', 2);
isempty(specs2)
nonlincols = strncmpi('nonlin',specs.Properties.VariableNames,6);
nonlinvars = table2cell(specs(:, nonlincols));
nonlinear = cell(size(nonlinvars,1),1);
for i = 1:size(nonlinvars,1)
    nonlinear{i} = strjoin(nonlinvars(i,:));
end
prefs = struct2cell(table2struct([nonlinear, specs(:, ~nonlincols)]));
