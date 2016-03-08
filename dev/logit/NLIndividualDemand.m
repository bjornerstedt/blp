classdef NLIndividualDemand < NLDemand
% NLIndividualDemand calculates individual choices    
    methods
        function q = getDemand(obj, p)
            q = zeros(size(obj.dummarket, 1), obj.settings.nind);
            for t = 1:size(obj.dummarket, 2)
                obj.initSimulation(t);
                sh = obj.shares(p(obj.dummarket(:, t)), t);
                cutoff1 = tril(ones(length(sh)),-1) * sh ;
                cutoff2 = tril(ones(length(sh))) * sh;
                qi = q(obj.dummarket(:, t), :);
                for i = 1:obj.settings.nind
                    r=rand;
                    qi(:,i) = cutoff1 < r & r < cutoff2;
                end
                q(obj.dummarket(:, t),:) = qi;
            end            
        end
        
        function obj = NLIndividualDemand(varargin)
            obj = obj@NLDemand(varargin{:});
            obj.settings.setParameters({'nind'});
            
            obj.settings.nind = 100;            
            obj.results.estimateDescription = 'Individual Nested Logit Demand';             
        end
    end
    
end

