classdef MixedLogitDemandMarket  < matlab.mixin.Copyable 
    % Random Coefficient demand market class
    %   Market level calculations for MixedLogitDemand.
    %   $Id:  $
    
    properties
        selection 
        config
        nind
        iweight
        v
        vx
        x2
        expmu
        s
        p
        d
        nonlinprice 
        alpha
        rc_sigma
    end
    
    methods         
        function [s, indshare] = shares(obj, p)
            if obj.config.ces
                p = log(p);
            end          
            if any(obj.nonlinprice) && any( p ~= obj.p)
                obj.p = p;
                obj.x2(:, obj.nonlinprice) = p;
                % Can update just nonlinprice
                for k = 1:size( obj.x2, 2)
                    obj.vx(:,:,k) = bsxfun(@times, obj.x2(:,k), obj.v(:,:,k));
                end
            end
            % obj.d set in init() or by user in simulating data:
            edeltapred = exp(obj.d - obj.alpha * p); 
            [s, indshare] = obj.sharecalc(edeltapred);
        end
       
        % One can test whether substituting the rhs for all instances of
        % obj.vx(:,:,k) affects performance. 
        function init(obj)
            for k = 1:size(obj.x2, 2)
                obj.vx(:,:,k) = bsxfun(@times, obj.x2(:,k), obj.v);
            end
            obj.nlpart(obj.rc_sigma);
        end
        
        function der = deltaJacobian(obj, rc_sigma, edelta)
            [~, si] = obj.sharecalc(edelta(obj.selection));
            wsi =  bsxfun(@times, si , obj.iweight'); % Row means
            dssigma = zeros(size(si, 1), length(rc_sigma));
            for k = 1:length(rc_sigma)
                svx = wsi .* obj.vx(:,:,k); % Can be done outside loop
                dssigma(:, k) = sum(svx, 2) - si * sum(svx)';
            end
            dsdelta = diag(sum(wsi, 2)) - wsi * si';
            % Note that dividing by obj.settings.nind cancels out.
            der = -dsdelta\dssigma;
        end
        
        function newedel = finddelta(obj, rc_sigma, edelta, tolerance)
            i = 0;
            maxdev = 100;
            obj.nlpart( rc_sigma );
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
        
        function nlpart(obj, rc_sigma)
            [n,m,K] = size(obj.vx);
            mu = zeros(n, m);
            for k = 1:K
                mu = mu + obj.vx(:,:,k) * rc_sigma(k);
            end
            obj.expmu = exp(mu);
        end
        
    end
end

