classdef NestedLogitDemand < DemandEstimate
    % Estimation and simulation code for nested logit.
    %   $Id: NestedLogitDemand.m 118 2015-05-26 10:34:48Z d3687-mb $
    
    properties
        alpha
        sigma1 = 0
        sigma2 = 0
    end
    properties % (SetAccess = protected )
        sigma
        nest1
        nest2
        sel % demand selection, used in selecting a market
        d % Utility delta
        xi %
        G % Group membership, each column a group with 1 if member
        GH % Subgroup membership
        GG % Group membership 
        GHGH % Subgroup membership 
        GN % Binary  of subgroup membership in group
    end
    
    methods
        function select(obj, selection)
            obj.sel = copy(obj); % Selection is set with select() method
            obj.sel.initSimulation(selection);
        end
        
        function createPeriods(obj)
            obj.period = cell(obj.nmkt,1);
            periodCopy = copy(obj); 
            periodCopy.T = [];
            for t = 1:max(obj.marketid)
                selection = obj.marketid == t;
                newperiod = copy(periodCopy); 
                newperiod.T = obj.T(selection, :);
                newperiod.initSimulation();
                obj.period{t} = newperiod;
            end            
        end
              
        function initSimulation(obj, varargin)
            if nargin > 1
                selection = varargin{1};
                obj.T( ~selection , :) = [];
            end
            obj.init();
            price = obj.getPrice();
            priceName = obj.getPriceName();
            shares = obj.generateShares(); 
            if isempty(obj.nestlist) 
                estimate = obj.results.estimate{priceName, 1}';
                obj.d = shares.ls - price * estimate';
            elseif length(obj.nestlist) == 1
                estimate = obj.results.estimate{[priceName, {'lsjg'}], 1}';
                obj.d = shares.ls - [price, shares.lsjg] * estimate';
            elseif length(obj.nestlist) == 2
                estimate = obj.results.estimate{[priceName, {'lsjh', 'lshg'}], 1}';
                obj.d = shares.ls - [price, shares.lsjh, shares.lshg] * estimate';
            end
            
            if ismatrix(estimate)
                obj.alpha = estimate(1);
                if length(estimate) > 1
                    obj.sigma = estimate(2:end);
                end
            end
			if  ~isempty(obj.nestlist)  
                obj.nest1 = obj.T{: , obj.nestlist(1)};
                gr = grp2idx(obj.nest1);
                obj.G = dummyvar(gr)';
				obj.GG = obj.G' * obj.G;
				if length(obj.nestlist) == 2 
                    obj.nest2 = obj.T{: , obj.nestlist(2)};
                    [~,~, g] = unique(obj.T(:, strsplit(strtrim(obj.var.nests))), 'rows');
                    obj.GH=dummyvar(g)';
   					obj.GHGH = obj.GH'*obj.GH;
					obj.GN =( (obj.G*obj.GH') > 0);
                    obj.GN = obj.G*obj.GH';
					obj.GN(obj.GN >0)=1;
				end
			end	
        end

        function price = getPrice(obj)
            if obj.settings.ces
                price = log(obj.T{:,obj.var.price});
            else
                price = obj.T{:,obj.var.price};
            end
        end

		function s = shares(obj, P)
            if obj.settings.ces
                P = log(P);
            end
			delta = obj.d + obj.alpha .* P; % d is delta without price effect
			if length(obj.nestlist) == 2 
                sigma1 = obj.sigma(1);
                sigma2 = obj.sigma(2);
				ev=exp(delta ./ (1-sigma1));
				ighs = (obj.GH*ev);
				igh = obj.GH' * ighs;
				evv = ighs .^ ((1-sigma1)/(1-sigma2));
				igs = (obj.GN * evv);
				ig = obj.G' * igs;
				it = sum(igs .^(1-sigma2));
				s = ev .* (igh .^ ((sigma2-sigma1)/(1-sigma2)) ) .* ...
                    (ig .^ (-sigma2)) / (1+it);
			elseif length(obj.nestlist) == 1 
                sigma1 = obj.sigma(1);
				ev=exp(delta ./ (1-sigma1));
				igs = (obj.G * ev);
				ig = obj.G' * igs;
				it = sum(igs .^(1-sigma1));
				s = ev .* (ig .^ (-sigma1)) / (1+it);
			else
				ev=exp(delta );
				it = sum(ev);
				s = ev / (1+it);
			end 
		end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
		% Used in Market.findcosts()
        function S = actualDemand(obj)
            S = obj.s;	
        end

		% Calculate quantities and market shares
        function tab = quantity(obj, P)
            tableCols = {'Product', 'Price', 'Quantity', 'MarketSh', 'Share'};
            s = obj.shares(P);
            if obj.settings.ces
                q = s .* obj.ms ./ P;
            else
                q = s .* obj.ms;
            end
            m = q / sum(q); 
            tab = table(obj.panelid, P, q, m, s);
            tab.Properties.VariableNames = tableCols;
        end
                
        function sj = shareJacobian(obj, P)
            if isempty(P)  
                S = obj.s;
                P = obj.p;
            else
                S = obj.shares(P);
            end
			other_effect =  - S*S';
			if length(obj.nestlist) == 0 	
				own_effect = diag(S);
				sj = obj.alpha *( own_effect + other_effect );  
            else
                sigma1 = obj.sigma(1);
            
				Sg = obj.GG * S;
				own_effect = 1/(1 - sigma1)*diag(S);
				if length(obj.nestlist) == 2 
                    sigma2 = obj.sigma(2);
					gr_effect = -sigma2/(1 - sigma2)*obj.GG .* ((S ./ Sg)*S');
					Sgh = obj.GHGH * S;
					gr_effect = gr_effect - (1/(1-sigma1) - 1/(1-sigma2))* ...
                        obj.GHGH .* ((S ./ Sgh)*S');
				else
					gr_effect = -sigma1/(1 - sigma1)*obj.GG .* ((S ./ Sg)*S');
				end 
				sj = obj.alpha*( own_effect + gr_effect + other_effect );  
			end
            if obj.settings.ces
                sj = (sj - diag(S) ) * diag(1./P) ;
            end
        end
        
        function s = sumstats(obj, j, E)
            e = j .* E;
            e(e==0)=[];
            s = [mean(e) std(e) min(e) max(e)];
        end
        
        function [elas, varargout] = elasticities(obj, P)
            s = obj.shares(P);
            n = length(s);
            D = obj.shareJacobian(P)';
            E = diag(P) * D * diag( 1 ./ s );
            elas = [obj.sumstats(eye(n), E)];
			if length(obj.nestlist) == 0 	
                elas = [elas; obj.sumstats(1 - eye(n), E)]
                rowtit = {'e_ii', 'e_ij'};
            elseif length(obj.nestlist) == 1 	
                elas = [elas; ...
                    obj.sumstats(obj.GG - eye(n), E); ...
                    obj.sumstats(1 - obj.GG, E)];
                rowtit = {'e_ii', 'e_ij', 'e_ik'};
            elseif length(obj.nestlist) == 2	
                elas = [elas; ...
                    obj.sumstats(obj.GHGH - eye(n), E); ...
                    obj.sumstats(obj.GG - obj.GHGH , E); ...
                    obj.sumstats(1 - obj.GG, E)];
                rowtit = {'e_ii', 'e_ij', 'e_ik', 'e_il'};
            end
            elas = array2table(elas);
            elas.Properties.VariableNames = {'Mean', 'Std', 'Min', 'Max'};
            elas.Properties.RowNames = rowtit;
            if nargout > 0
                varargout{1} = E;
            end
        end
        
        function [elas] = groupElasticities(obj, P,  group)
            group = obj.T{:, group};
            if iscategorical(group)
                [names,~,group] = unique(group);
                names = matlab.lang.makeValidName(cellstr(char(names)));
            else
                names = [];
            end
            A = dummyvar(group)';
            s = obj.shares(P);
            D = obj.shareJacobian(P)';
            elas = A*diag(P)* D * A' * diag( 1 ./(A*s) ) ;
            elas = array2table(elas);
            if ~isempty(names)
                elas.Properties.RowNames = names;
                elas.Properties.VariableNames = names;
            end
        end
        
        function obj = NestedLogitDemand(varargin)
            if nargin > 0
                obj.T = varargin{1};
            end
            obj.settings.estimateMethod = 'gls';
            obj.config.estimateDescription = 'Nested Logit Demand'; 
        end
    end
end

