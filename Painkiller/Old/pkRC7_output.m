function  pkRC_output( input, testdata, job, test  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fvals = [];
for n = 1: size(job,2)
    if isempty(job{n})
        continue;
    end
    if n == 1
        results = table2array(job{n}.estimate2);
    else
        results(:,:,n) = table2array(job{n}.estimate2);
    end
    coefsnew = job{n}.estimate2(:, 1);
    coefsnew.Properties.VariableNames =  {sprintf('Sim%d', n)};
    stddevnew = job{n}.estimate2(:, 2);
    stddevnew.Properties.VariableNames =  {sprintf('Sim%d', n)};
    mergernew = job{n}.merger(:, 4);
    mergernew.Properties.VariableNames =  {sprintf('Sim%d', n)};
    costsnew = job{n}.merger(:, 1);
    costsnew.Properties.VariableNames =  {sprintf('Sim%d', n)};
    fvals = [fvals job{n}.fval];
    if n == 1
        coefs = coefsnew;
        stddev = stddevnew;
        costs = costsnew;
        merger = mergernew;
    else
        coefs = [coefs, coefsnew];
        stddev = [stddev, stddevnew];
        merger = [merger, mergernew];
        costs = [costs, costsnew];
    end
end
resAv = mean(results,3);
resStd = std(results,0,3);
fvaltable=array2table(fvals);
fvaltable.Properties.RowNames = {'fval1', 'fval2','optimalIt','cond', 'sdreg'}';
fvaltable.Properties.VariableNames = coefs.Properties.VariableNames ;
coefsff = union(coefs,fvaltable,'Rows','stable')
if ispc
    writetable(coefsff,input.fn,'WriteRowNames',true,'Sheet', test  );
    loc = sprintf('A%d', 28);
    writetable(stddev,input.fn,'WriteRowNames',true, 'Sheet', test, 'Range',loc );
    loc = sprintf('A%d', 51);
    writetable(merger,input.fn,'WriteRowNames',true,'Sheet', test, 'Range',loc );
    loc = sprintf('A%d', 61);
    writetable(costs,input.fn,'WriteRowNames',true, 'Sheet', test, 'Range',loc );
else
    writetable(coefsff,'results_coef.csv','WriteRowNames',true );
    writetable(stddev,'results_stddev.csv','WriteRowNames',true );
end

end

