clear
x = [-2, .3; 0, .4; 2, .3; ];
nonlinear = {{'test'}, {'testemp', 'empirical', x}, {'testuni', 'uniform'}};

draws = Draws(...
    'DrawMethod', 'quadrature',...
    'Individuals', 5,...
    'Accuracy', 3 ...
    );
draws.parse( nonlinear)

draws.create( );

% draws.draws
% 
% draws.weights
