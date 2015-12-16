classdef Sampling
    %SAMPLING Various sampling methods
    %   General methods with static invocation
    
    methods(Static)
        function w = draw(method, N, draws, varargin)
            if nargin > 2 
                rs = varargin{1};
            else
                rs = [];
            end
            switch lower(method)
                case 'hypercube'
                    w = Sampling.mlhs(draws, N, rs);
                case 'halton'
                    w = Sampling.halton(N, draws);
                case 'halton2'
                    w = Sampling.drawhalton(N, draws);
                case 'random'
                    if isempty(rs)
                        w = randn(draws, N);
                    else
                        w = rs.randn(draws, N);
                    end
                otherwise
                    error('Incorrect settings.drawmethod specification')
            end
        end
                    
        function draws = mlhs( N, D, varargin)
            if nargin > 2 
                rs = varargin{1};
            else
                rs = [];
            end
            if ~isempty(rs)
                ordereddraws = zeros(N,D);
                shuffleddraws = zeros(N,D);
                for i = 1:D
                    ordereddraws(:,i) = ((1:N)-1)/N + rs.rand()/N;
                    shuffle = rs.randperm(N);
                    shuffleddraws(:,i) = ordereddraws(shuffle,i);
                end
                draws = norminv(shuffleddraws);
            else
                ordereddraws = zeros(N,D);
                shuffleddraws = zeros(N,D);
                for i = 1:D
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

