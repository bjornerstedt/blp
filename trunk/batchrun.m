function [job, diary]  = batchrun(func, data, testdata, varargin)
%BATCHRUN Execute a set of tests
%   To use, implement a script and a function.
%   The script: 
%       1) initializes a cell array of structs: testdata
%       2) initializes a struct of common data: data
%       3) initializes a cell array of random streams if 'randomstream' is
%          specified
%       4) invokes func( data, testdata, repetition) or if 'randomstream'
%          is spec func( data, testdata, repetition, randstream)
%       5) processes the result, also a cell array of structs
%
%   The function func receives the following parameters: 
%        1) a struct with common data 
%        1) a struct with test specific data 
%        1) the repetition number 
%        1) a random stream rs specific for the repetition. Use rs.rand() to
%           make the randomization in the reptition independent from
%           others.
%   the current test number. It returns all data as struct, put in a TxN cell
%   array by batchrun, where T is the length of testdata and N is the 
%   repetitions. It should not have any random number creation, with 
%   individual and common data put in structures testdata and data 
%   respectively. For example nind has to be created and included in the 
%   data struct. The indata from testdata and data are combined by 
%   batchrun to a struct, sent to func. The returned cell array diary
%   contains the output of the execution of func in the case of parallel
%   processing.
%   Invoking batchrun(... , 'process',@savefunc) invokes this function
%   after each test. Parameters passed are the test struct, the cell array of
%   results and the test number. Can be used for tabulation or other
%   processing.
%   $Id: batchrun.m 111 2015-04-30 11:18:10Z d3687-mb $

p = inputParser;
p.addRequired('func', @(x)isa(x, 'function_handle'));
p.addRequired('data',@isstruct);
p.addRequired('testdata',@isstruct);
p.addOptional('repetitions',1,@isscalar);
p.addParameter('process',[], @(x)isa(x, 'function_handle'));
p.addParameter('parallel',false, @islogical);
p.addParameter('select',[], @ismatrix);
p.addParameter('randomstream', 0, @isscalar);
p.parse(func, data, testdata, varargin{:});

if p.Results.randomstream ~= 0
    randstream = RandStream.create('mrg32k3a', ...
        'NumStreams', p.Results.repetitions,...
        'CellOutput',true,'Seed', p.Results.randomstream);
end
% Set randstream{i} as stream in each repetition or if not set, presuppose that
% RandStream.setGlobalStream(s2) has been invoked with some seed.

N = p.Results.repetitions;
T = length(testdata);
job = cell(T,N);
if ~isempty(p.Results.select)
    % A single test has been selected with 'select' option
    if T == 1 && length(p.Results.select) == 1
        N=p.Results.select;
    elseif length(p.Results.select) == 2
        T=p.Results.select(1);
        N=p.Results.select(2);
    else
        error('Selection must be a scalar when 1 test exists and two elements when there are several tests.');
    end
    test_T = feval(func, p.Results.data, p.Results.testdata(T), 0,[]);
    job = {feval(func, p.Results.data, test_T, N, randstream{N})};
    if ~isempty(p.Results.process)
        feval(p.Results.process, data, test_T, job, 1);
    end
    diary = '';
    return
end
if T > 1
    fprintf('Evaluating %d tests, repeating %d times\n', T, N);
else
    fprintf('Executing single test %d times\n', N );
end
if ~p.Results.parallel || (T == 1 && N == 1)
    for t = 1:T
        test_t = feval(func, p.Results.data, p.Results.testdata(t), 0, []);
        for n = 1:N
            if p.Results.randomstream ~= 0
                job{t,n}  = feval(func, p.Results.data, test_t, n, randstream{n});
            else
                job{t,n}  = feval(func, p.Results.data, test_t, n);
            end
        end
        if ~isempty(p.Results.process)
            feval(p.Results.process, data, test_t, job(t,:), t);
        end
    end
    diary = {};
else
    diary = cell(T,1);
    pool = gcp();
    timeout = 30*60; % Max time per eval 30 min
    for t = 1:T
        test_t = feval(func, p.Results.data, p.Results.testdata(t), 0, []);
        
        for n = 1:N
            if p.Results.randomstream ~= 0
                f(n) = parfeval(pool, @dojob,1,func, p.Results.data, test_t, n,...
                    randstream{n}); 
            else
                f(n) = parfeval(pool, @dojob,1,func, p.Results.data, test_t, n); 
            end
        end
        parfor_progress(N);
        for n = 1:N
            [completedIdx,value] = fetchNext(f, timeout);
            if ~isempty(value)
                job{t, completedIdx} = value;
            end
            parfor_progress;
        end
        parfor_progress(0);
        diary{t} = '';
        for n = 1:N
            disp(f(n).Diary)
            diary{t} = [diary{t}, f(n).Diary];
        end
        if ~isempty(p.Results.process)
            feval(p.Results.process, data, test_t, job(t,:), t);
        end
    end
end
end

function job = dojob(func, varargin)
    try
        job = feval(func, varargin{:});
    catch err
        warning('batchrun caught error', err.identifier);
        job = [];
    end
end
