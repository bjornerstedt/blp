% **************************************************************
% pkRC7_run.m runs a single estimation with variations in simulation
% Output is rendered by pkRC_output.m
% **************************************************************
clear
tic
input.repetitions = 10;
input.fn = 'test.xls';
input.marketdraws = true;

nonlinear = {'paracetamol fizzy branded'; 'paracetamol fizzy branded constant'};
ces = {true; true};
commonDraws  = {true; true};
testdata = table2struct(table(nonlinear, ces, commonDraws));
[job,diary] = batchrun(@pkRC7, input, testdata, input.repetitions, ...
    'save', @pkRC_output, 'Parallel', true);

rt = toc
rt/3600