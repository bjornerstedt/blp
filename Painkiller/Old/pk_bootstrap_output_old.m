% load 'mergedata'
headings = {'Firm', 'Firm2', 'Product' , 'Costs', 'NewCosts' ,'NC' , 'Merge', 'MergeNC', 'MergeNCD', 'Bootrep' };
evalcols = {'NC' , 'Merge', 'MergeNC', 'MergeNCD', 'Costs', 'NewCosts'};
mergedata.Properties.VariableNames = headings;
fname = table(mergedata.Firm);
[Firm,~,idx]= unique(mergedata.Firm);
mergedata.firm = idx;
[Firm2,~,idx]= unique(mergedata.Firm2);
mergedata.firmsubst = idx;
mergedata.Firm = [];
mergedata.Firm2 = [];
mergedata.Product = [];
if bootreps >= 100;
    grvars = {'Bootrep','firm'};
    fmean =  varfun(@mean, mergedata, 'GroupingVariables', grvars);
    fmean.GroupCount = [];
    fmean.mean_firmsubst = [];
    fmean.Properties.VariableNames = [grvars, evalcols];
    pos=1
    for n = 1:length(evalcols)
        prop = evalcols{n};
        nd = sortrows(fmean, {'firm', prop});
        prods = size(nd,1)/bootreps;
        res = reshape(nd{:,prop}, bootreps, prods);
        perc = [1 , 5, 50, 95, 99]' ;
        perclev = perc * bootreps/100;
        bootlevels = array2table([perc, res(perclev,:)]);
        bootlevels.Properties.VariableNames = [{'Percentage'}; cellstr(char(Firm))];
        disp ' '
        disp(['Confidence levels for: ', prop]);
        disp(bootlevels)
        if saveXLS
            xlswrite(fn, {['Confidence levels for: ', prop]}, 'intervals', ...
                sprintf('A%d', pos) );
            writetable(bootlevels, fn, 'Sheet', 'Intervals', 'Range', ...
                sprintf('A%d', pos+2) );
        end
        pos = pos + 10;
    end
end

fname = {table(Firm),table(Firm2)};
func = {@mean, @std, @min, @max};
label = {'Mean', 'Standard deviation', 'Min', 'Max'};
xlspos = {'A%d', 'J%d'};
firmvar = {'firm', 'firmsubst'};
pos = 1;
if saveXLS
    xlswrite(fn, {'Results for estimate'}, 'results', 'A1' );
    writetable(standardMerger, fn, 'sheet', 'results','WriteRowNames', true, 'range', 'A3');
    pos = pos + 12;
end
for n = 1:length(func)
    for r = 1:length(firmvar)
        tab =  varfun(func{n}, mergedata, 'GroupingVariables', firmvar{r});
        tab = [fname{r}, tab(:,3:8)];
        tab.Properties.VariableNames = [{'Firm'}, evalcols];
        if saveXLS
            xlswrite(fn, label(n), 'results', sprintf(xlspos{r}, pos) );
            writetable(tab, fn, 'sheet', 'results', 'range', sprintf(xlspos{r}, pos+2));
        end
    end
    pos = pos + 12;
end
if saveXLS
    writetable(mergedata, fn, 'sheet', 'bootdata');
    writetable(bootdraws, fn, 'sheet', 'draws');    
    writetable(alldata, fn, 'sheet', 'AllData');    
    if false
    writetable(demand.results.settings, fn,'WriteRowNames',true, ...
        'Sheet', 'settings');
    writetable(demand.results.var ,filename,'WriteRowNames',true, ...
        'Sheet', 'settings', 'Range', 'E1' );    
    end
end
