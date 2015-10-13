% **************************************************************
% run_pkRC.m runs a single estimation with variations in simulation
% Output is rendered by pkRC_output.m
% **************************************************************
clear
tic
input.repetitions = 2;
input.fn = 'resultsquad.xls';

% nonlinear = cell(2,1);
% nonlinear(:) = {'paracetamol constant'};
% ces = {1; 1};
% commonDraws = {true; true};
% marketdraws = {false; true};
% nind  = [500; 500];
nonlinear = cell(1,1);
nonlinear(:) = {'paracetamol fizzy branded constant'};
ces = { 1};
commonDraws = { true};
marketdraws = {true};
nind  = [ 500];

testdata = table2struct(table(nonlinear, ces, commonDraws, ...
    marketdraws, nind));

% [job,diary] = batchrun(@pkRC, input, testdata, input.repetitions, ...
%     'save', @pkRC_output, 'Parallel', false,'select', 2);
[job,diary] = batchrun(@pkRC, input, testdata, input.repetitions, ...
    'save', @pkRC_output, 'Parallel', false);

rt = toc
rt/3600