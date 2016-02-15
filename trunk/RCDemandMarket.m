classdef RCDemandMarket  < matlab.mixin.Copyable 
    % Market calculations in Mixed Logit
    %   Used by RCDemand 
    
    properties
        selection 
        config
        iweight
        v
        vx
        x2
        expmu
        s
        d
        nonlinprice 
        alpha
        sigma
    end
    
    methods         
        function obj = RCDemandMarket(demand)
            % Copy config as a struct to avoid deep copy problems with
            % copy() function. Note that one could keep the common config
            % object, since all instances have same config. 
            obj.config = demand.config.getProperties();
            obj.alpha = demand.alpha;
            obj.sigma = demand.sigma;
            obj.config.ces = demand.settings.ces;
            obj.nonlinprice = strcmp(demand.var.price, demand.nonlinparams);    
        end
        
        function [s, indshare] = shares(obj, p)
            if obj.config.ces
                p = log(p);
            end          
            if any(obj.nonlinprice) 
                k = obj.nonlinprice;
                obj.vx(:,:,k) = bsxfun(@times, p, obj.v(k,:));
            end
            % obj.d set in init() or by user in simulating data:
            edeltapred = exp(obj.d - obj.alpha * p); 
            [s, indshare] = obj.sharecalc(edeltapred);
        end
       
        function sh = shareJacobian(obj, P)
            [S, si] = obj.shares(P);
            if any(obj.nonlinprice) 
                sigma_p = obj.sigma(obj.nonlinprice);
                % obj.sim.v is the same for all products
                dUi = obj.iweight .* (- obj.alpha + ...
                    sigma_p * obj.v( obj.nonlinprice, :))';
            else
                dUi = - obj.alpha * obj.iweight;
            end
            sh = zeros(size(si,1));
            for i = 1:size(si,2)
                shi = si(:, i);
                sh = sh + dUi(i) * (diag(shi) - shi * shi');
            end
            if obj.config.ces
                sh = (sh - diag(S) ) * diag(1 ./ P) ;
            end            
        end
        
        % One can test whether substituting the rhs for all instances of
        % obj.vx(:,:,k) affects performance. 
        function init(obj)
            for k = 1:size(obj.x2, 2)
                obj.vx(:,:,k) = bsxfun(@times, obj.x2(:,k), obj.v(k,:));
            end
            obj.nlpart(obj.sigma);
        end
        
        function der = deltaJacobian(obj, sigma, edelta)
            % Note that sigma is an unecessary parameter. The following
            % line does not affect the estimation outcome.
            obj.nlpart(sigma);
            [~, si] = obj.sharecalc(edelta(obj.selection));
            wsi =  bsxfun(@times, si , obj.iweight'); % Row means
            dssigma = zeros(size(si, 1), length(sigma));
            for k = 1:length(sigma)
                svx = wsi .* obj.vx(:,:,k); % Can be done outside loop
                dssigma(:, k) = sum(svx, 2) - si * sum(svx)';
            end
            dsdelta = diag(sum(wsi, 2)) - wsi * si';
            % Note that dividing by obj.settings.nind cancels out.
            der = -dsdelta\dssigma;
        end
        
        function newedel = findDelta(obj, sigma, edelta, tolerance)
            i = 0;
            maxdev = 100;
            obj.nlpart( sigma );
            while maxdev > tolerance && i < obj.config.fpmaxit
                eg = bsxfun(@times, edelta, obj.expmu );
                % Code with intermediate steps
                % denom = 1 ./ ( 1 + sum(eg));
                % share = bsxfun(@times, eg, denom) * obj.iweight; 
                % newedel = edelta .* obj.s ./ share;
                newedel = edelta .* obj.s ./ ( ...
                    bsxfun(@times, eg, 1 ./ ( 1 + sum(eg)) ) ...
                    * obj.iweight);
                maxdev = max(abs(newedel - edelta));       
%                     dev = newedel - edeltasel;   % Is this motivated?
%                     maxdev = dev'*dev;
                edelta = newedel;
                i = i + 1;
            end
        end
        
        function [sh, indsh] = sharecalc(obj, edelta)
            eg = bsxfun(@times, edelta, obj.expmu );
            denom = 1 ./ ( 1 + sum(eg));
            indsh = bsxfun(@times, eg , denom);
            sh = indsh * obj.iweight; % Row means
        end
        
        function nlpart(obj, sigma)
            [n,m,K] = size(obj.vx);
            mu = zeros(n, m);
            for k = 1:K
                mu = mu + obj.vx(:,:,k) * sigma(k);
            end
            obj.expmu = exp(mu);
        end
        
    end
end

