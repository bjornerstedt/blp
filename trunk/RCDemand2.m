classdef RCDemand2 < RCDemand   
    % Add Hessian to RCDemand
    properties
    end
    
    methods
        function [f, g, h] = objective(obj, sigma)
            %This function defines the objective over which to minimize
            obj.edelta =  obj.findDelta(sigma);
            del = log(obj.edelta);
            
            if max(isnan(del)) == 1
                f = 1e+10;
                if nargout > 1
                    g = 1e+10;
                end
            else
                [~, xi] = obj.lpart(del);
                if ~isempty(obj.Z)
                    f = xi' * obj.ZWZ * xi;
                    if nargout == 2
                        obj.deltaJac = obj.deltaJacobian(sigma, obj.edelta);
                        g = 2* obj.deltaJac'* obj.ZWZ * xi;
                    elseif nargout == 3
                        [g, h] = obj.objectiveHessian(xi, del, sigma);
                    end
                else
                    if false
                        f = xi' * xi;
                        if nargout == 2
                            obj.deltaJac = obj.deltaJacobian(sigma, obj.edelta);
                            g = 2*obj.deltaJac' * xi;
                        elseif nargout == 3
                            [g, h] = obj.objectiveHessian(xi, del, sigma);
                        end
                    else
                        f = xi' * xi;
                        if nargout == 2
                            obj.deltaJac = obj.deltaJacobian(sigma, obj.edelta);
%                             g = 2*obj.deltaJac'* obj.X * xiX';
                        g = 2*obj.deltaJac' * xi;
                        elseif nargout == 3
                            [g, h] = obj.objectiveHessian(xi, del, sigma);
                        end
                    end
                end
            end
        end
 
        function [ bet, xi ] = lpart(obj, del )
%             if strcmpi(obj.settings.paneltype, 'fe')
%                 davt = accumarray(obj.panelid, del, [], @mean);
%                 davt = davt(obj.panelid, :);
%                 del = del - davt;
%             end
            if ~isempty(obj.Z)
                bet = obj.inv_x1ZWZx1 * (obj.X' * obj.ZWZ * del);
            else
                bet = obj.inv_x1ZWZx1 * (obj.X' * del);
            end
%             bet = obj.estimationMatrix * del;

            xi = del - obj.X * bet;
        end
        
        function [j, h] = objectiveHessian(obj, xi, delta, sigma)
            % Create vector expression from Hessian
            % $-2ZWZ^{T}\xi(\theta)$
            if ~isempty(obj.Z)
                hessVector = 2 * obj.ZWZ * xi;
            else
                hessVector = 2 * xi;
            end
            hess = zeros(length(sigma),length(sigma));
            jac = zeros(size(delta));
            j2 = zeros(size(sigma));
            for t = 1:max(obj.marketid)
                [jp, hp] = obj.period{t}.objectiveHessian( delta(obj.dummarket(:,t)), sigma, ...
                    hessVector(obj.dummarket(:,t)) );
                hess = hess + hp;
                jac(obj.dummarket(:,t)) = jp; 
                j2 = j2 + jp'* hessVector(obj.dummarket(:,t));
            end
            j = jac' * hessVector;
            if ~isempty(obj.Z)
                h = 2*jac' * obj.ZWZ * jac - hess;
            else
                h = 2 * jac' * jac - hess;
            end
        end
 
        function sigma = minimize(obj, varargin)
            %        extraoptions = {'OutputFcn', @iterOutput}; % Show output func
            if nargin > 1
                extraoptions = varargin{1};
            else 
                extraoptions = {};
            end
            if obj.config.hessian
                extraoptions = [extraoptions, {'Hessian', 'on'}];
            end
            if isempty(obj.sigma)
                error('sigma is not set, probably because init has not been invoked');
            end
            obj.oldsigma = zeros(size(obj.sigma));
            options = optimoptions(@fminunc, ...
                'MaxIter', obj.settings.maxiter, ...
                'TolX', obj.config.tolerance, 'TolFun', obj.config.tolerance, ...
                'Algorithm','trust-region' ,...
                'GradObj','on', ...
                extraoptions{:});
            func = @(x)obj.objective(x);
            [sigma,fval,exitflag,output,grad,hessian] = fminunc(func, obj.sigma, options);
        end
        
        function obj = RCDemand2(varargin)
            obj = obj@RCDemand(varargin{:});
        end
%             obj.var.setParameters({'nonlinear','nonlinearlogs','nonlineartriangular'});
%             obj.settings.setParameters({'optimalIV','drawmethod',... 
%                 'marketdraws','nind','quaddraws','maxiter','fptolerance1','fptolerance2'});
%             obj.config = SettingsClass({'tolerance','fpmaxit', ...
%                 'restartMaxIterations','restartFval', 'test', ...
%                 'quadrature','randstream','guessdelta','quietly'});
%         end      

        function initPeriods(obj)
            md = RCDemandMarket(obj);
            md.iweight = obj.draws.weights;
            obj.period = cell(max(obj.marketid), 1);
            for t = 1:max(obj.marketid)
                newmarket = copy(md);
                newmarket.selection = logical(obj.dummarket(:, t));
                newmarket.x2 = obj.x2(newmarket.selection, :);
                if ~isempty(obj.share)
                    newmarket.s = obj.share.s(newmarket.selection, 1);
                end
                if ~isempty(obj.d)
                    newmarket.d = obj.d(newmarket.selection, 1);
                end
                newmarket.v = [];
                if obj.settings.marketdraws
                    newmarket.v = obj.draws.draws(:,:,t)'; 
                else
                    newmarket.v = obj.draws.draws'; 
                end
                newmarket.init();
                obj.period{t} = newmarket;
            end
        end
        
        function newedelta = findDelta(obj, sigma)
            if isempty(obj.edelta)
                % Create starting values for findDelta
                if isempty(obj.share.s)
                    obj.share.s = sum(bsxfun(@times, obj.dummarket, ...
                        1 ./ sum(obj.dummarket) ), 2) .* 0.5;
                    obj.edelta = obj.share.s .* 2;
                else
                    obj.edelta = obj.share.s ./ obj.share.s0;
                end
            end
            if max(abs(sigma - obj.oldsigma)) < 0.01
                tolerance = obj.config.fptolerance2;
                closeFlag = 0;
            else
                tolerance = obj.config.fptolerance1;
                closeFlag = 1;
            end
            
            if obj.config.guessdelta && ~isempty(obj.deltaJac)
                newdelta = log(obj.edelta) + ...
                    obj.deltaJac*(sigma - obj.oldsigma); 
                edelta = exp(newdelta);
            else
                edelta =  obj.edelta;
            end

            sel = logical(obj.dummarket);
            newedelta = zeros(size(edelta));
            for t = 1:max(obj.marketid)
                newedelta(sel(:,t)) = obj.period{t}.findDelta(sigma, ...
                    edelta(sel(:,t)), tolerance);
            end
            % Update oldsigma and edelta only in first stage and if
            % successful
            if closeFlag == 1 && max(isnan(newedelta)) < 1;
                obj.oldsigma = sigma;
            end
        end       
    end
end

