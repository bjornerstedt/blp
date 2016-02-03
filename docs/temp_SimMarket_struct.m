%% SimMarket linear parameter structure
% beta does not include alpha:
beta = [ 1, 0];

% x is the expected values, const not included. The value for p only
% matters if market is calculated, not simulated. 
x = [5, 0];
x_vcv = [1, 1];
% Intercept, should be called gamma
c = 4;
% Corresponds to epsilon_sigma
c_vcv = 1;
gamma = 0;

% Individual and product level shocks:
epsilon_sigma = .1;
sigma_xi = .1;

endog_sigma = 0.1; % Degree of corr between x and epsilon
prob_prod = .8;    % Prob of product existing in market    

epsilon = endog_sigma + sigma_xi
eta = N( 0, c_vkv)

x = N(model.x, model.x_vcv)
d = [x(2:end), 1] * beta' + epsilon

w = N(0, 1)
c = modell.c + model.gamma * w + eta