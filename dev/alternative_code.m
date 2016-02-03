        function newedelta = findDeltaExact(obj, rc_sigma, edelta)
            selection = logical(obj.dummarket);
            obj.nlpart(rc_sigma); % Creates obj.expmu
            results = cell(max(obj.marketid),1);
            options = optimoptions(@fsolve,'Display','off');
            newedelta = edelta;
            for t = 1:max(obj.marketid)
                edeltasel = edelta(selection(:,t),1);
                expmu = obj.expmu(selection(:,t), :);
                sorig = obj.s(selection(:,t),1);
                edeltasol  = fsolve(...
                    @(x)(sorig - obj.marketShare(x, expmu)), ...
                    edeltasel, options);
                newedelta(selection(:,t)) = edeltasol;
            end
        end
    
        function share = marketShare(obj, edelta, expmu)
            %Computes the sum of the individual shares
            eg = bsxfun(@times, edelta, expmu );
            denom = 1 ./ ( 1 + sum(eg));
            if obj.quadrature
                share = bsxfun(@times, eg, denom) * obj.quadw; % Row means
            else
                share = mean(bsxfun(@times, eg, denom), 2); % Row means
            end
        end
  
        function fval = minimize2(obj, varargin)
            extraoptions = {'OutputFcn', @iterOutput};
            extraoptions = {};
            if nargin > 1
                extraoptions = [extraoptions, varargin{1}];
            else
                %extraoptions = {};
            end
            obj.historyx = [];
            obj.historyfval = [];
            obj.oldsigma = zeros(size(obj.x2,2),1);

            %This function estimates the random coefficient model
            if obj.fminsearch
                options = optimset( ...
                    'TolX', obj.tolerance, 'TolFun', obj.tolerance, extraoptions{:});
                %,'MaxFunEvals',1000,'MaxIter',1000);
                [rc_sigma,fval,exitflag,output] = fminsearch( ...
                    @objective, obj.rc_sigma, options);
            else
                options = optimoptions(@fminunc, ...
                    'TolX', obj.tolerance, 'TolFun', obj.tolerance, ...
                    'GradObj','on', extraoptions{:});
                %,'MaxFunEvals',1000,'MaxIter',1000);
                [rc_sigma,fval,exitflag, output] = fminunc( ...
                    @objective, obj.rc_sigma, options);
            end
            obj.results.fval = fval;
            obj.rc_sigma = rc_sigma;
            
            %%%%% Sub functions
            function [f, g] = objective(rc_sigma)
                %This function defines the objective over which to minimize
                if max(abs(rc_sigma - obj.oldsigma)) < 0.01
                    tolerance = obj.fptolerance2;
                    flag = 0;
                    obj.startDeltaZero = true;
                else
                    tolerance = obj.fptolerance1;
                    flag = 1;
                end
                %%%%%%%%%%%%%%%%%
                i = 0;
                maxContractionAttempts = 2;
                while(i < maxContractionAttempts) 
                    if obj.startDeltaZero
                        obj.edelta = ones(size(obj.edelta));
                    end
                    if obj.exactDelta
                        edel = obj.findDeltaExact(rc_sigma, obj.edelta);
                    else
                    if obj.guessdelta && ~isempty(obj.deltaJac)
                        newdelta = log(obj.edelta) + ...
                            obj.deltaJac*(rc_sigma - obj.oldsigma);
                        edel = obj.finddelta(rc_sigma, exp(newdelta), tolerance);
                    else
                        edel = obj.finddelta(rc_sigma, obj.edelta, tolerance);
                    end
                    end

                    if any(isnan(edel))
                        obj.startDeltaZero = true
                        i = i+1;
                    else
                        i = maxContractionAttempts;
                    end
                end
                %%%%%%%%%%%%%%%%%
                % Update oldsigma and edelta only in first stage and if
                % successful
                if any(isnan(edel))
                    warning('delta from contraction is NaN');
                end
                if flag == 1 && max(isnan(edel)) < 1;
                   obj.oldsigma = rc_sigma;
                end   
                
                obj.edelta = edel;
                del = log(edel);
                
                if max(isnan(del)) == 1
                    f = 1e+10;
                    if nargout > 1
                        g = 1e+10;
                    end
                else
                    [bet, xi] = obj.lpart(del);
                    f = xi' * obj.ZWZ * xi;
                    if nargout > 1
                        obj.deltaJac = obj.deltaJacobian(rc_sigma, obj.edelta);
                        g = 2* obj.deltaJac'* obj.ZWZ * xi;
                    end
                end
            end
            
            function stop = iterOutput(x, optimValues, state)
                stop = false;
                switch state
                    case 'iter'
                        fprintf('Iterations: %d / %d  \t Sigma: %2.4f\n', ...
                            obj.fp_invoc, obj.fp_iter, x(1));
                        obj.historyx = [obj.historyx ; x'];
                        obj.fp_iter = 0;
                        obj.fp_invoc = 0;
                end
            end
        end
        


% Objective function with repeated attempts to find fixedpoint
            function [f, g] = objective(rc_sigma)
                %This function defines the objective over which to minimize
                if max(abs(rc_sigma - obj.oldsigma)) < 0.01
                    tolerance = obj.fptolerance2;
                    flag = 0;
                else
                    tolerance = obj.fptolerance1;
                    flag = 1;
                end
                %%%%%%%%%%%%%%%%%
                i = 0;
                maxContractionAttempts = 5;
                while(i < maxContractionAttempts) 
                if obj.exactDelta
                    edel = obj.findDeltaExact(rc_sigma, obj.edelta);
                else
                if obj.guessdelta && ~isempty(obj.deltaJac)
                    newdelta = log(obj.edelta) + ...
                        obj.deltaJac*(rc_sigma - obj.oldsigma);
                    edel = obj.finddelta(rc_sigma, exp(newdelta), tolerance);
                else
                    edel = obj.finddelta(rc_sigma, obj.edelta, tolerance);
                end
                end
                
                if any(isnan(edel))
                    obj.edelta = exp( log(obj.edelta) + randn(size(obj.edelta)));
                    i = i+1;
                else
                    i = maxContractionAttempts;
                end
                end
                %%%%%%%%%%%%%%%%%
                % Update oldsigma and edelta only in first stage and if
                % successful
                if any(isnan(edel))
                    warning('delta from contraction is NaN');
                end
                if flag == 1 && max(isnan(edel)) < 1;
                   obj.oldsigma = rc_sigma;
                end   
                
                obj.edelta = edel;
                del = log(edel);
                
                if max(isnan(del)) == 1
                    f = 1e+10;
                    if nargout > 1
                        g = 1e+10;
                    end
                else
                    [bet, xi] = obj.lpart(del);
                    f = xi' * obj.ZWZ * xi;
                    if nargout > 1
                        obj.deltaJac = obj.deltaJacobian(rc_sigma, obj.edelta);
                        g = 2* obj.deltaJac'* obj.ZWZ * xi;
                    end
                end
            end
            
            
