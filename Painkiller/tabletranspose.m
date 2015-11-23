function tabout = tabletranspose( tabin )
% Transpose of numeric table
    tabout = array2table(table2array(tabin)');
    tabout.Properties.VariableNames = tabin.Properties.RowNames;
    tabout.Properties.RowNames = tabin.Properties.VariableNames;
end

