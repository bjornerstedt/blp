function  pkRC_output( input, testdata, job, test  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

est = length(job{1}.demand);
coefs = cell(est,1);
stddev  = cell(est,1);
properties  = cell(est,1);
merger = [];
costs = [];
rc_sigma0 = [];
for n = 1: size(job,2)
    if isempty(job{n})
        continue;
    end
    
    colh = {sprintf('Sim%d', n)};
    for i = 1:est
        
        column = job{n}.demand(i).results.estimate(:, 1);
        column.Properties.VariableNames = colh;
        coefs{i} = [coefs{i}, column];
        column = job{n}.demand(i).results.estimate(:, 2);
        column.Properties.VariableNames = colh;
        stddev{i} = [stddev{i}, column];
        
        column = job{n}.demand(i).results.properties;
        column.Properties.VariableNames = colh;
        properties{i} = [properties{i}, column];
    end
    column = job{n}.merger(:, 4);
    column.Properties.VariableNames = colh;
    merger = [merger, column];
    
    column = job{n}.merger(:, 1);
    column.Properties.VariableNames = colh;
    costs = [costs, column];
    rc_sigma0 = [rc_sigma0 job{n}.demand(1).results.rc_sigma0];
end
rc_sigma0 = array2table(rc_sigma0);
rc_names = coefs{1}.Properties.RowNames(end-size(rc_sigma0,1)+1:end);
rc_sigma0.Properties.VariableNames = coefs{1}.Properties.VariableNames;
rc_sigma0.Properties.RowNames = rc_names;

if isempty(coefs)
    warning('No simulations to display');
    return;
end
display 'Coefficients'
disp(coefs{est});

if input.save
if ispc 
    pos = 3;
    result1 = job{1}.demand(1).results;
    filename = [input.fn,'.xls'];
    xlswrite(filename, {'Settings'}, test, 'A1'); 
    writetable(result1.settings, filename,'WriteRowNames',true, ...
        'Sheet', test, 'Range',sprintf('A%d', pos) );
    writetable(result1.var ,filename,'WriteRowNames',true, ...
        'Sheet', test, 'Range',sprintf('E%d', pos) );
    pos = pos + max(size(result1.settings,1), ...
        size(result1.var,1));
    
    pos = outtable('Coefficients', coefs{est}, filename, test, pos  );
    pos = outtable('Properties First stage', properties{1}, filename, test, pos  );
    if length(properties) == 2
        pos = outtable('Properties', properties{2}, filename, test, pos  );
    end
    pos = outtable('Standard Deviations', stddev{est}, filename, test, pos  );
    pos = outtable('Price Change', merger, filename, test, pos  );
    pos = outtable('Costs', costs, filename, test, pos  );
    pos = outtable('Starting Points', rc_sigma0, filename, test, pos  );

else
    writetable(coefs{est},'results_coef.csv','WriteRowNames',true );
    writetable(stddev{est},'results_stddev.csv','WriteRowNames',true );
end
end

if true
    for i = 1:est
        coefs{i} = tabletranspose( coefs{i} );
        stddev{i} = tabletranspose( stddev{i} );
        properties{i} = tabletranspose( properties{i} );
    end
    rc_sigma0 = tabletranspose( rc_sigma0);
    price_ch = tabletranspose( merger );
    costs = tabletranspose( costs);
    demand = job{1}.demand(est);
    filename = [input.fn,'.mat'];
    save(filename, 'demand', 'coefs', 'stddev', 'rc_sigma0', ...
        'price_ch', 'costs');
end

end

function pos = outtable(  heading, table, fn, test, pos  )
    spacing = 4;
    xlswrite(fn, {heading}, test, sprintf('A%d', pos+2)); 
    pos = pos + spacing;
    writetable(table , fn, 'WriteRowNames',true, ...
        'Sheet', test, 'Range', sprintf('A%d', pos) );
    pos = pos + size(table, 1);
end

function tabout = tabletranspose( tabin )
% Transpose of numeric table
    tabout = array2table(table2array(tabin)');
    tabout.Properties.VariableNames = tabin.Properties.RowNames;
    tabout.Properties.RowNames = tabin.Properties.VariableNames;
end

