% *********************************************
% pkRC_run.m is the main script to run pkRC.m:
% * several tests
% * several repetitions
% Tests are run in parallel or sequentially
% Output is rendered by pkRC_output.m
% *********************************************
clear
tic
input.repetitions = 1
input.fn = 'test1.xls';
input.ces = true;
input.commonDraws = true;
input.marketdraws = true;

if input.ces
    price = 'Ptablets ';
else
    price = 'Ptablets_Real ';
end
commonRC = 'constant ' % RC in all treatments
nonlinear = {'paracetamol '; ...
    'paracetamol fizzy '};
% ['paracetamol fizzy branded ']; ...
% ['paracetamol fizzy ', price ]; ...
% ['paracetamol ',commonRC]; ...
% ['paracetamol fizzy ',commonRC];...
% ['paracetamol fizzy branded ',commonRC]; ...
% ['paracetamol fizzy ', price ,commonRC]};

ces = {true; true};
commonDraws  = {true; true};
testdata = table2struct(table(nonlinear, ces, commonDraws));
[job,diary] = batchrun(@pkRC, input, testdata, input.repetitions, ...
    'save', @pkRC_output, 'Parallel', false);

rt = toc
rt/3600