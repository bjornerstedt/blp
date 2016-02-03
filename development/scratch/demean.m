function f  = demean( T, vars, indexvars )
%DEMEAN Subtract average value of vars
%   The index variables are the variables averages are taken over.
    [uniqueRows,~,rowIdx]=unique(T{:,indexvars},'rows');
    T = T( :, [indexvars vars]);
    B = varfun(@mean, T,'GroupingVariables', indexvars);
    B =  B{rowIdx, length(indexvars)+2:end}
    B = T{:, length(indexvars)+1:end} - B ;
    B = array2table(B);
    B.Properties.VariableNames = vars;
    f = B;
end
