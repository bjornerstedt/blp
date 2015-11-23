% demand = MixedLogitDemand();
% xlswrite( 'test.xls', [fieldnames(demand.settings),struct2cell(demand.settings)])

var = fieldnames(demand.var);
value = struct2cell(demand.var);
vartable = table(var,value)
setting = fieldnames(demand.settings);
value = struct2cell(demand.settings);
settingstable = table(setting,value)
writetable( settingstable, 'test.xls', 'Sheet', 'Settings')
writetable( vartable, 'test.xls', 'Sheet', 'Settings', 'Range', 'E1')
