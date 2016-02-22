clear all
SimMarket.randDraws();

%% Test 1: empirical draws
distlist = {{'test'}, {'testuni', 'uniform'}, {'testuni', 'triangular'}};
dr1 = Draws('DrawMethod', 'hypercube', 'individuals', 1e6);
dr1.create( distlist);

for j = 1:length(distlist)
    SimMarket.testSame(max(mean(dr1.draws)), 0)
    SimMarket.testSame(max(var(dr1.draws) - 1), 0)
end

display ' **** Test 1 passed ****'

%% Test 2: empirical draws
dr1 = Draws('DrawMethod', 'quadrature', 'Accuracy', 3);
dr1.create( {{'test'}, {'testuni', 'uniform'}});

dr2 = Draws('DrawMethod', 'quadrature', 'Accuracy', 3);
[dm, wt] = nwspgr('KPN', 1, 3);
dr2.create( {{'testemp', 'empirical', [dm, wt]}, {'testuni', 'uniform'}});
assert(all(all(dr1.draws - dr2.draws == 0)))
assert(all(all(dr1.weights - dr2.weights == 0)))
display ' **** Test 1 passed ****'

%% Test 3: moments of random draws
distname = {'normal', 'uniform', 'empirical', 'triangular', ...
    'logistic', 'lognormal'};
drawmethod = {'hypercube', 'halton', 'random'};

for i = 1:length(drawmethod)
    results = zeros(length(distname), 2);
    
    for j = 1:length(distname)
        if j== 3
            results(j, 1) = 0;
            results(j, 2) = 1;
            continue
        end
        nonlinear = {{'test1', distname{j} }};
        
        draws = Draws(...
            'DrawMethod', drawmethod{i},...
            'Individuals', 1000000 ...
            );
        draws.create(nonlinear);
        
        results(j, 1) = mean(draws.draws);
        results(j, 2) = var(draws.draws);
    end
    % Max difference for distributions:
    results - repmat([0, 1], size(results,1) ,1 )
    max(max(abs(results - repmat([0, 1], size(results,1) ,1 ))))
    
    results = array2table(results);
    results.Properties.VariableNames = {'mean', 'std'};
    results.Properties.RowNames = distname;
    
    disp(['Drawmethod: ', drawmethod{i}])
    disp(results)
end

%% Test 4: moments of quadrature draws
distname = {'normal', 'uniform'};

results = cell(length(distname), 2);

for j = 1:length(distname)
    if j== 3
        results{j, 1} = 0;
        results{j, 2} = 1;
        continue
    end
    nonlinear = {{'test1', distname{j} }};
    
    draws = Draws(...
        'DrawMethod', 'quadrature',...
        'Accuracy', 7 ...
        );
    
    draws.create( nonlinear);
    
    results{j, 1} = mean(draws.draws);
    % Variance calculated as weighted summation:
    results{j, 2} = sqrt(sum(draws.draws .^ 2 .* draws.weights));
end

results = cell2table(results);
results.Properties.VariableNames = {'mean', 'std'};
results.Properties.RowNames = distname;

disp('Drawmethod: quadrature')
disp(results)

%% Test 5: Market draws
distname = {'normal', 'uniform', 'empirical', 'triangular', ...
    'logistic', 'lognormal'};
drawmethod = {'hypercube', 'halton', 'random'};

draws = Draws(...
    'DrawMethod', 'hypercube',...
    'Markets', 2, ...
    'Individuals', 1000 ...
    );
draws.create(  'test1');
dr = reshape(draws.draws, 1, []);
mean(dr)
std(dr)