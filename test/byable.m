% Create table
A = [1 2 3 6;1 2 2 5;1 4 1 4;2 2 1 8;2 2 3 1];
T = array2table(A);
T.Properties.VariableNames = {'a', 'b', 'c', 'd'};

index = {'a', 'b'};

B = varfun(@mean, T,'GroupingVariables', index)

% Remove index and group count
B(:,index) = [];
B.GroupCount = [];
B.Properties.RowNames = {}; % Not necessary whan converting to array
B

[uniqueRows,~,rowIdx]=unique(T{:, index},'rows');
B =  B{rowIdx, :}
T(: , index) = [];
table2array(T) - B % Cannot subtract tables
break

%% USING accumarray instead

%# find unique rows and their corresponding indices in A
[uniqueRows,~,rowIdx]=unique(A(:,1:2),'rows');

%# for each group of unique rows, sum the values of the third column of A
count = accumarray(rowIdx,A(:,3),[],@length);

subtotal =accumarray(rowIdx,A(:,3))
[A subtotal(rowIdx,:) ]

% A collapse is 
[uniqueRows,subtotal,count]

av = accumarray(rowIdx,A(:,3:end),[],@mean);
av(rowIdx,:)