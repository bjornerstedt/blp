
filename = [fn,'.xlsx'];
bootreps = max(alldata{:,'Bootrep'});
changevar = {'Price' , 'WeightedPrice' , 'MarketSh' , 'Costs' };
changevarfun = {@mean, @sum, @sum, @mean};
resultcat = {'Original', 'Coll', 'Merge', 'MergeNC', 'MergeColl', 'MergeNCcoll'};
% resultcat = {'Original', 'Merge', 'MergeNC'};
% resultcat(2) = { 'Coll', 'MergeColl', 'MergeNCcoll'};
commonvars = {'Firm', 'Firm2', 'Product', 'Bootrep', 'Substance'};
grvar = {'Substance', 'Firm2'};
origShares = alldata{alldata.Result==resultcat{1} & alldata.Bootrep==1, 'MarketSh'};
origShares = repmat(origShares,size(alldata,1)/size(origShares,1),1);
alldata.WeightedPrice = alldata.Price .* origShares;
evalcols = resultcat(3:end);
excelcol = {'A','G','M','S'};
for c = 1:length(changevar)
    pos=1;
    for g = 1:length(grvar)
        grtemp = unique(alldata(:,grvar{g}));
        grcats = cellstr(char(grtemp{:,grvar{g}}));
        grvars = {'Bootrep', 'Result', grvar{g}};
        summary =  varfun(changevarfun{c}, alldata, 'GroupingVariables', grvars, ...
            'InputVariables', changevar{c});
        summary.GroupCount = [];        
        % Compare percentage changes over resultcat with resultcat{1}
        for rc = 3:length(resultcat)
            results = summary(summary.Result == resultcat{1}, grvars);
            v1 = summary{summary.Result == resultcat{1}, end};
            for r = 3:length(resultcat)
                v2 = summary{summary.Result == resultcat{r}, end};
                name = resultcat{r};
                if c ~= 3
                    results{:, name} = (v2 - v1) ./ v1;
                else
                    results{:, name} = (v2 - v1);
                end
            end
        end
        if bootreps < 100;
            results(:,{'Bootrep', 'Result', grvar{g}}) = [];
            results.Properties.RowNames = grcats;
            disp(results);
            writetable(results, filename, 'Sheet', changevar{c}, 'Range', ...
                sprintf('A%d', pos), 'WriteRowNames', true );
            pos = pos + size(results, 1) +4;
        else
            for n = 1:length(evalcols)
                meanval =  varfun(@mean, results, 'GroupingVariables', grvar{g}, ...
                    'InputVariables', evalcols{n});
                
                nd = sortrows(results, {grvar{g}, evalcols{n}});
                prods = size(nd,1)/bootreps;
                res = reshape(nd{:,evalcols{n}}, bootreps, prods);
                perc = [50, 10, 90]' ;
                perclev = round(perc * bootreps/100);
                bootlevels = [meanval(:,end),array2table(res(perclev,:)')];
                bootlevels.Properties.VariableNames = ...
                    {'Mean', 'Median', 'CI_Lower', 'CI_Upper'};
                bootlevels.Properties.RowNames = grcats;
                disp ' '
                outtext = {grvar{g}, changevar{c}};
                if display
                    disp(outtext);
                    disp(bootlevels)
                end
                xlswrite(filename, outtext, evalcols{n}, ...
                    sprintf('%s%d', excelcol{c}, pos) );
                writetable(bootlevels, filename, 'Sheet', evalcols{n}, 'Range', ...
                    sprintf('%s%d', excelcol{c}, pos+1), 'WriteRowNames', true );
            end
            pos = pos + size(bootlevels, 1) +4;
        end
    end
end

