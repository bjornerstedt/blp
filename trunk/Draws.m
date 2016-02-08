classdef Draws < matlab.mixin.Copyable
    %SAMPLING Various sampling methods
    %   General methods with static invocation
    properties
        settings
        
        names
        draws
        weights
    end
    methods
        function draw(obj, K)
            nind = obj.settings.individuals;
            obj.weights = ones(nind, 1) / nind ;
            drawcount = nind * obj.settings.markets;
            rs = obj.settings.randstream;
            switch lower(obj.settings.drawmethod)
                 case 'hypercube'
                    dm = Draws.mlhs(drawcount, K, rs);
                case 'halton'
                    dm = Draws.halton(K, drawcount);
                case 'halton2'
                    dm = Draws.drawhalton(K, drawcount);
                case 'random'
                    if isempty(rs)
                        dm = randn(drawcount, K);
                    else
                        dm = rs.randn(drawcount, K);
                    end
                otherwise
                    error('Incorrect settings.drawmethod specification')
            end
            dm = dm';
            if obj.settings.markets > 1
                obj.draws = [];
                for t = 1:obj.settings.markets
                    obj.draws(:, :, t) = dm(: , (t - 1)*nind + (1:nind));
                end
            else
                obj.draws = dm;
            end
        end
        
        function quadrature(obj, K, acc)
            [X, obj.weights] = nwspgr('KPN', K, acc);
            obj.draws = X';
        end
        
        function obj = Draws(varargin)
            args = inputParser;
            args.addParameter('drawmethod', 'hypercube', @ischar);
            args.addParameter('accuracy', 7, @isnumeric);
            args.addParameter('individuals', 300, @isnumeric);
            args.addParameter('markets', 1, @isnumeric);
            args.addParameter('randstream', []);
            args.parse(varargin{:});
            obj.settings = args.Results;
        end
        
    end
    methods(Static)
                            
        function draws = mlhs( K, D, varargin)
            if nargin > 2 
                rs = varargin{1};
            else
                rs = [];
            end
            if ~isempty(rs)
                ordereddraws = zeros(K,D);
                shuffleddraws = zeros(K,D);
                for i = 1:D
                    ordereddraws(:,i) = ((1:K)-1)/K + rs.rand()/K;
                    shuffle = rs.randperm(K);
                    shuffleddraws(:,i) = ordereddraws(shuffle,i);
                end
                draws = norminv(shuffleddraws);
            else
                ordereddraws = zeros(K,D);
                shuffleddraws = zeros(K,D);
                for i = 1:D
                    ordereddraws(:,i) = ((1:K)-1)/K + rand()/K;
                    shuffle = randperm(K);
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
        
        function w = halton( vars, draws)
            truncdist = 10;
            hs = haltonset( vars, 'Skip', 10);
            w = net(hs, draws);
            w = norminv(w);
            w(w > truncdist) = truncdist;
            w(w < -truncdist) = -truncdist;
        end
        
        % Replication of Mathias draws
        function hm = drawhalton(vars, draws)
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

