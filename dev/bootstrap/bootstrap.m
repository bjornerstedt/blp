% bootstrap returns a bootstrap draw of a table, grouped by indeces.
function data2 = bootstrap(data, indeces)
    % create vector of unique indexes
    [~, g, c] = unique(data{:,indeces});
    g2 = datasample(g,length(g));
    data2 = table();
    for i = 1:length(g)
        newdata = data(c==g2(i),:);
        newdata{:,indeces} = i;
        data2 = [data2; newdata]; 
    end
end
