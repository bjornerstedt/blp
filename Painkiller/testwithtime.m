%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RC bootstrap with cost and demand calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% Primary Config parameters

RC = true    % Set RC to choose between RC and NL
ces = true   % Set ces to choose between ces and unit demand
select = false  % Choose between multiple RC starting points or selection  

% Other parameters

estimate = true
display = true;
saveXLS = true;

input.repetitions = 20;
input.save = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nonlinear = cell(2,1);
nonlinear(:) = {'paracetamol fizzy branded constant'};
nind  = [500; 500];
cesdemand = {true; true};
withtime = {true; true};
newinstruments = {false; true};

testdata = table2struct(table(nonlinear, cesdemand, nind, withtime, newinstruments));
if select
    input.nocons = true;
    input.fn = 'test1';
    input.packDemand = false;
    [job,diary] = batchrun(@pkRCwithtime, input, testdata, input.repetitions, ...
        'Parallel', false, 'select', selectRun, ...
        'randomstream', 99);
    demand = job{1}.demand;
else
    input.nocons = false;
    input.fn = 'pk RC CES with time.xlsx';
    input.packDemand = true;
    [job,diary] = batchrun(@pkRCwithtime, input, testdata, input.repetitions, ...
        'process', @pkRC_output, 'Parallel', true, 'randomstream', 99);
end
