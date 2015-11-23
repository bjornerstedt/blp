% *********************************************
% Main script to run either 
% * several repetitions of a test
% * all tests
% Tests are run in parallel
% *********************************************
clear
tic
if true % merges job when there are many starting points in the same test
    saveResults = false
    for test = 1:2
    repetitions = 2
    global commondraws; % For saving individual draws
    global repcount; % For keeping track of the number of repetitions
    repcount = 1;
    % [job,diary] = dobatch(@pkML1, test, repetitions);
     [job,diary] = dobatchseq(@pkML1, test, repetitions);
   
    fvals = [];
    for n = 1: size(job,1)
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
    fvaltable.Properties.RowNames = {'fval1', 'fval2','optimalIt','cond','isreal', 'sdreg'}';
    fvaltable.Properties.VariableNames = coefs.Properties.VariableNames ;
    coefsff = union(coefs,fvaltable,'Rows','stable')
    if saveResults
        fn = 'results_ces1.xls';
        %  fn = 'results_unit2.xls';
        if ispc
            writetable(coefsff,fn,'WriteRowNames',true,'Sheet', test  );
            loc = sprintf('A%d', 28);
            writetable(stddev,fn,'WriteRowNames',true, 'Sheet', test, 'Range',loc );
            loc = sprintf('A%d', 51);
            writetable(merger,fn,'WriteRowNames',true,'Sheet', test, 'Range',loc );
            loc = sprintf('A%d', 61);
            writetable(costs,fn,'WriteRowNames',true, 'Sheet', test, 'Range',loc );
        else
            writetable(coefsff,'results_coef.csv','WriteRowNames',true );
            writetable(stddev,'results_stddev.csv','WriteRowNames',true );
        end
    end
    end
else
    % merges job when there are many tests
    [job,diary] = dobatch(@pkML);
    for n = 1: size(job,1)
        if n == 1
            coefs = job{n}.estimate2(:, 1);
            %   coefs = job{n}.merger(:, 1);
            coefs.Properties.VariableNames = {sprintf('Sim%d', n)};
            coefs.keys=coefs.Properties.RowNames;
        else
            restabnew = job{n}.estimate2(:, 1);
            %   restabnew = job{n}.merger(:, 1);
            restabnew.Properties.VariableNames =  {sprintf('Sim%d', n)};
            restabnew.keys=restabnew.Properties.RowNames;
            coefs = outerjoin(coefs, restabnew,'MergeKeys',true);
        end
    end
    coefs.Properties.RowNames = coefs.keys;
    coefs.keys=[];
    disp(coefs)
    
    if saveResults
    % On pc do write excel
        writetable(coefs,'results.csv','WriteRowNames',true );
    end
end
rt = toc
rt/3600