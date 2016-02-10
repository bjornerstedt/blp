clear

%% Test 1: empirical draws
x = [-1, .5; 1, .5 ];
nonlinear = {{'test'}, {'testemp', 'empirical', x}, {'testuni', 'uniform'}};
draws = Draws(...
    'DrawMethod', 'hypercube',...
    'Individuals', 500,...
    'Accuracy', 3 ...
    );
draws.parse( nonlinear)

draws.create( );
mean(draws.draws)
var(draws.draws)

%% Test 2: moments of random draws
distname = {'normal', 'uniform', 'empirical', 'triangular', ...
    'logistic', 'lognormal'};
drawmethod = {'hypercube', 'halton', 'random'};

for i = 1:length(drawmethod)
    results = cell(length(distname), 2);

for j = 1:length(distname)
if j== 3
    results{j, 1} = 0;
    results{j, 2} = 1;
    continue
end
nonlinear = {{'test1', distname{j} }};

draws = Draws(...
    'DrawMethod', drawmethod{i},...
    'Individuals', 5000 ...
    );
draws.parse( nonlinear);
draws.create( );

results{j, 1} = mean(draws.draws);
results{j, 2} = var(draws.draws);
end

results = cell2table(results);
results.Properties.VariableNames = {'mean', 'std'};
results.Properties.RowNames = distname;

disp(['Drawmethod: ', drawmethod{i}])
disp(results)
end

%% Test 3: moments of quadrature draws
draws = Draws(...
    'DrawMethod', 'hypercube',...
    'Individuals', 50000,...
    'Accuracy', 3 ...
    );
draws.parse( nonlinear)

draws.create( );
mean(draws.draws)
var(draws.draws)