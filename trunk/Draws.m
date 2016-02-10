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
        function [dm, wt] = generate(obj, distribution, K)
            nind = obj.settings.individuals;
            wt = ones(nind, 1) / nind ;
            drawcount = nind * obj.settings.markets;
            rs = obj.settings.randstream;
            switch lower(obj.settings.drawmethod)
                case 'hypercube'
                    dm = Draws.mlhs(drawcount, K, rs);
                case 'halton'
                    dm = Draws.halton(drawcount, K);
                case 'halton2'
                    dm = Draws.drawhalton(drawcount, K);
                case 'random'
                    % Note that these are normal rather than uniform draws
                    if isempty(rs)
                        dm = randn(drawcount, K);
                    else
                        dm = rs.randn(drawcount, K);
                    end
                otherwise
                    error('Incorrect settings.drawmethod specification')
            end
            if obj.settings.markets > 1
                obj.draws = [];
                for t = 1:obj.settings.markets
                    obj.draws(:, :, t) = dm(:, (t - 1)*nind + (1:nind));
                end
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
            [X, wt] = nwspgr(distmeth, K, obj.settings.accuracy);
            if distribution == 2
                % Standardize uniform distribution
            end
            dm = X;
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
        
        function create2(obj)
            % Simple creation process if obj.var.nonlinear is a string
            K = length(obj.names);
            if strcmpi(obj.settings.drawmethod, 'quadrature')
                [obj.draws, obj.weights] = obj.quadrature( 1, K);
            else
                [obj.draws, obj.weights] = obj.generate( 1, K);
            end
            obj.draws = obj.draws';
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
                [dl{1, 1}, dl{1, 2}] = obj.generate( 2, length(vars) );
            end
            % Put in a local function and invoke for each quad/rand
            obj.draws = dl{1, 1};
            obj.weights = dl{1, 2};
            for i = 2:size(dl, 1)
                [obj.draws, obj.weights] = Draws.gridCombine(...
                    obj.draws, obj.weights, dl{i, 1}, dl{i, 2});
            end
            % Market draws comes here
            obj.draws = obj.draws';
        end
        
        function names = parse(obj, spec)
            distnames = {'normal', 'uniform', 'empirical'};
            if isstr(spec)
                names = strsplit(strtrim(spec));
                obj.names = names;
            elseif iscell(spec)
                % Handle random draws differently: combine all but
                % empirical to uniform, and inverse map CDF on each based
                % on some kind of list.
                names = {};
                for i = 1:length(spec)
                    speci = spec{i};
                    namesi = strsplit(strtrim(speci{1}));
                    names =[ names, namesi];
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
                            
        function draws = mlhs( N, K, varargin)
            if nargin > 2 
                rs = varargin{1};
            else
                rs = [];
            end
            if ~isempty(rs)
                ordereddraws = zeros(N,K);
                shuffleddraws = zeros(N,K);
                for i = 1:K
                    ordereddraws(:,i) = ((1:N)-1)/N + rs.rand()/N;
                    shuffle = rs.randperm(N);
                    shuffleddraws(:,i) = ordereddraws(shuffle,i);
                end
                draws = norminv(shuffleddraws);
            else
                ordereddraws = zeros(N,K);
                shuffleddraws = zeros(N,K);
                for i = 1:K
                    ordereddraws(:,i) = ((1:N)-1)/N + rand()/N;
                    shuffle = randperm(N);
                    shuffleddraws(:,i) = ordereddraws(shuffle,i);
                end
                draws = norminv(shuffleddraws);
            end
        end
        
        function draws = triangular( unif)
            low = sqrt(6)*(sqrt(2*unif) - 1);
            high = sqrt(6)*(1 - sqrt(2*(1 - unif)) );
            low(unif>=0.5) = 0;
            high(unif<0.5) = 0;
            draws = high+low;
        end
        
        function w = halton( draws, vars)
            truncdist = 10;
            hs = haltonset( vars, 'Skip', 10);
            w = net(hs, draws);
            w = norminv(w);
            w(w > truncdist) = truncdist;
            w(w < -truncdist) = -truncdist;
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
        
        % Replication of Mathias draws
        function hm = drawhalton(draws, vars)
            prim = [109 11 101 11 13 103 41 43 47 79 97 73 107 ...
                89 3 7 13 17 19 23 29 31 37 53 59 61 67 71 83 113];
            % prim=[2 3 5];
            hm=[];
            for h = 1:vars
                hm1=halton(10+(draws), prim(h));
                hm1=norminv(hm1);
                %CENSOR normal halton draws to 10 and -10
                hm1=hm1.*(hm1 <= 10)+10.*(hm1>10);
                hm1=hm1.*(hm1 >= -10)-10.*(hm1<-10);
                help=hm1(11:size(hm1,1),1);
                hm=[hm reshape(help',draws,1)];
            end
        end
        
    end   
end

