classdef Draws < matlab.mixin.Copyable
% Draws - Various sampling methods
% The Draws class provides sampling methods for RCDemand. It is created as
% a separate class to separate this task from estimation, facilitating 
% creation of unit tests. 
    properties
        settings
        
        names
        spec
        draws
        weights
    end
    methods
        function [dm, wt] = generate(obj, varTypes, K)
            function d = normalDist(d)
                d = norminv(d);
                % Truncation previously only used by halton
                truncdist = 10;
                d(d > truncdist) = truncdist;
                d(d < -truncdist) = -truncdist;
            end
            function d = uniformDist(d)
                % Draws should have zero expectation and variance 1
                d = 2*sqrt(3) * (d - 1/2);
            end
            function d = triangularDist( d)
                % Symmetric triangular distribution
                low = sqrt(6) * ( sqrt(2 * d) - 1 );
                high = sqrt(6) * ( 1 - sqrt(2 * (1 - d)) );
                low(d >= 0.5) = 0;
                high(d < 0.5) = 0;
                d = high + low;
            end
            function d = logisticDist(d)
                d = -sqrt(3)/pi * log(1 ./ d - 1);
            end
            function d = lognormalDist(d)
                d = norminv(d);
                mu2 = (1 + sqrt(5))/2;
                sigma = sqrt(log(mu2));
                d =  exp(sigma * d) - sqrt(mu2);
            end
            invCDF = {@normalDist, @uniformDist, [], ... 
                @triangularDist, @logisticDist, @lognormalDist};
            % Transform uniform to distribution using inverse CDF
            function dm = invTransform(dm)
                if isempty(varTypes)
                    dm = normalDist(dm);
                else
                    for i = 1:length(varTypes)
                        func = invCDF{varTypes{i}};
                        dm(:, i) = func(dm(:, i));
                    end
                end
            end
            
            nind = obj.settings.individuals;
            wt = ones(nind, 1) / nind ;
            switch lower(obj.settings.drawmethod)
                case 'hypercube'
                    dm = obj.mlhs(nind, K);
                case 'halton'
                    dm = Draws.halton(nind, K);
                case 'random'
                    dm = obj.random(nind, K);
                otherwise
                    error('Incorrect settings.drawmethod specification')
            end
            % Market draws
            if obj.settings.markets > 1
                dn = [];
                dsh = rand(obj.settings.markets, size(dm, 2));
                for t = 1:obj.settings.markets
                    nv = bsxfun(@plus, dm, dsh(t, :));
                    dn(:, :, t) = invTransform(nv - floor(nv));
                end
                dm = dn;
            else
                dm = invTransform(dm);
            end
        end
                
        function [dm, wt] = quadrature(obj, distribution, K)
            switch distribution
                case 1
                    distmeth = 'KPN';
                case 2
                    distmeth = 'KPU';
                otherwise
                    error('Method not supported for quadrature.')
            end
            [dm, wt] = nwspgr(distmeth, K, obj.settings.accuracy);
            if distribution == 2
                % Standardize uniform distribution
                dm = 2*sqrt(3) * (dm - 1/2);
            end
        end
        
        function [dm, wt] = empiricalDraws(obj, drwt)
            if ~isnumeric(drwt) || any(size(drwt) <= 1)
                error('Empirical weights must be a numerical matrix')
            end
            dm = drwt(:, 1:end-1);
            wt = drwt(:, end);
            if abs(sum(wt) - 1) > 1e-6
                 error('Weights must be frequencies summing to 1')
            end 
        end
        
        function dm = random(obj, N, K)
            rs = obj.settings.randstream;
            if isempty(rs)
                dm = rand(N, K);
            else
                dm = rs.rand(N, K);
            end
        end
        
        function shuffleddraws = mlhs(obj, N, K)
            rs = obj.settings.randstream;
            if ~isempty(rs)
                ordereddraws = zeros(N,K);
                shuffleddraws = zeros(N,K);
                for i = 1:K
                    ordereddraws(:,i) = ((1:N)-1)/N + rs.rand()/N;
                    shuffle = rs.randperm(N);
                    shuffleddraws(:,i) = ordereddraws(shuffle,i);
                end
            else
                ordereddraws = zeros(N,K);
                shuffleddraws = zeros(N,K);
                for i = 1:K
                    ordereddraws(:,i) = ((1:N)-1)/N + rand()/N;
                    shuffle = randperm(N);
                    shuffleddraws(:,i) = ordereddraws(shuffle,i);
                end
            end
        end
        
        function create2(obj)
            % Simple creation process if obj.var.nonlinear is a string
            K = length(obj.names);
            if strcmpi(obj.settings.drawmethod, 'quadrature')
                [obj.draws, obj.weights] = obj.quadrature( 1, K);
            else
                [obj.draws, obj.weights] = obj.generate( [], K);
            end
        end
        
        function create(obj)
            if isempty(obj.spec)
                obj.create2()
                return
            end
            dl = {};
            % for each obj.spec create row in drawlist
            if strcmpi(obj.settings.drawmethod, 'quadrature')
                for i = 1:length(obj.spec)
                    speci = obj.spec{i};
                    if speci{2} == 3 % empirical draws
                        [dl{i, 1}, dl{i, 2}] = obj.empiricalDraws(speci{3});
                    else
                        [dl{i, 1}, dl{i, 2}] = obj.quadrature(speci{2}, ...
                            length(speci(1)));
                    end
                end
            else
                % Collect all draws but empirical, storing in a local list
                vars = {};
                varTypes = {};
                for i = 1:length(obj.spec)
                    speci = obj.spec{i};
                    if speci{2} == 3 % empirical draws
                        [dl{2, 1}, dl{2, 2}] = obj.empiricalDraws(speci{3});
                    else
                        vars = [vars, speci{1}];
                        varTypes = [varTypes, repmat(speci(2), size(speci(1)))];
                    end
                end
                [dl{1, 1}, dl{1, 2}] = obj.generate( varTypes, length(vars) );
            end
            % Put in a local function and invoke for each quad/rand
            obj.draws = dl{1, 1};
            obj.weights = dl{1, 2};
            for i = 2:size(dl, 1)
                [obj.draws, obj.weights] = Draws.gridCombine(...
                    obj.draws, obj.weights, dl{i, 1}, dl{i, 2});
            end
            % Market draws comes here
        end
        
        function names = parse(obj, spec)
            distnames = {'normal', 'uniform', 'empirical', 'triangular', ...
                'logistic', 'lognormal'};
            if ischar(spec) 
                names = strsplit(strtrim(spec));
                obj.names = names;
            elseif iscell(spec)
                % With one dist type convert to nested cell array if it is
                % not:
                if ischar(spec{1})
                    spec = {spec};
                end
                names = {};
                for i = 1:length(spec)
                    speci = spec{i};
                    speci{1} = strsplit(strtrim(speci{1}));
                    names =[ names, speci{1}];
                    if length(speci) == 1
                        speci{2} = 1;
                    else
                        speci{2} = find(strcmpi( speci{2}, distnames));
                    end
                    if speci{2} == 3 && length(speci) ~= 3
                        error('Empirical dist must specify values as third arg.')
                    end
                    obj.spec{i} = speci;
                end
            else
                error('Input to nonlinear draws must be a string or a cell array')
            end
            obj.names = names;
        end

        function obj = Draws(varargin)
            args = inputParser;
            args.addParameter('drawmethod', 'hypercube', @ischar);
            args.addParameter('accuracy', 7, @isnumeric);
            args.addParameter('individuals', 5, @isnumeric);
            args.addParameter('markets', 1, @isnumeric);
            args.addParameter('randstream', []);
            args.parse(varargin{:});
            obj.settings = args.Results;
        end
        
    end
    methods(Static)
                                          
        function w = halton( draws, vars)
            hs = haltonset( vars, 'Skip', 10);
            w = net(hs, draws);
        end

        function [X, W] = gridCombine(X1, W1, X2, W2)
            % gridCombine creates combined weighted draws from two sets of
            % draws. All combination of rows of X1 and X2 are generated.
            % Assumes that all matrices are column matrices.
            function X = matCombine(X1, X2)
                v1 = repmat(X1', size(X2, 1), 1 );
                v1 =  reshape(v1, size(X1, 2), [])';
                X2 = repmat(X2, size(X1, 1), 1 );
                X = [v1, X2];
            end
            X = matCombine(X1, X2);
            W = prod(matCombine(W1, W2), 2);
        end
        
        function [v, quadw] = initquadrature0(K, N)
            % This function is from revision 110 to generate and integrate
            % K different quadratures with weighting
            [X,W] = nwspgr('GQN', 1, N);
            v = [];
            % Quadrature with accuracy N and K vars => N^K ind:
            quadw = ones( 1, N^K);
            for k = 1:K
                % X' is row vector of ind
                vt = repmat(X', N^(K-k), N^(k-1) );
                v(:,k) =  reshape(vt, 1, N^K);

                wt = repmat(W', N^(K-k), N^(k-1) );
                quadw = quadw .* reshape(wt, 1, N^K);
            end
            quadw = quadw'; % Facilitates weighted sum
        end
        
    end   
end

