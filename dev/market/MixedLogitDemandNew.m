classdef MixedLogitDemandNew < MixedLogitDemand    
    methods      
        function [s, si] = shares(obj, p)
            [s, si] = obj.period{obj.sim.market}.shares( p );
        end
       
        function [sh, indsh] = sharesAll(obj, p)
            sel = logical(obj.dummarket);
            sh = zeros(size(p));
            indsh = zeros(length(p), obj.settings.nind);
            for t = 1:max(obj.marketid)
                [s, si] = obj.period{t}.shares( p(sel(:,t)) );
                sh(sel(:,t)) = s;
                indsh(sel(:,t)) = si;
            end
        end
       
        function initSimulation(obj, market)       
            % General init, move to estimate? Should NestedLogit have
            % similar code?
            obj.sim.market = market; % HACK to get market based routines to work
            obj.alpha = - obj.beta(1);
            if isempty(obj.d)
                % Create starting values for finddelta
                if isempty(obj.share.s)
                    obj.share.s = sum(bsxfun(@times, obj.dummarket, ...
                        1 ./ sum(obj.dummarket) ), 2) .* 0.5;
                    ed = obj.share.s .* 2;
                else
                    ed = obj.share.s ./ obj.share.s0;
                end
                edelta = obj.finddelta(obj.rc_sigma, ed,...
                    obj.settings.fptolerance1);
                obj.d = log(edelta) + obj.alpha*obj.X(:, 1);
                for t = 1:max(obj.marketid)
                    obj.period{t}.d = obj.d(obj.dummarket(:, t));
                end
            end
        end
        
        % Invoked from objective
        function der = deltaJacobian(obj, rc_sigma, edelta)
            der = zeros(size(obj.x2));
            for t = 1:max(obj.marketid)
                index = obj.marketid == t;
                der(index, :) = obj.period{t}.deltaJacobian(rc_sigma, edelta);
            end
        end
        
        function init(obj, varargin)
            init@MixedLogitDemand(obj, varargin{:});
            obj.initPeriods();
        end
        
        function initPeriods(obj)
            md = MixedLogitDemandMarket();
            % Copy config as a struct to avoid deep copy problems with
            % copy() function:
            md.config = obj.config.getProperties();
            md.nind = obj.settings.nind;
            md.alpha = obj.alpha;
            md.rc_sigma = obj.rc_sigma;
            md.config.ces = obj.settings.ces;
            md.iweight = obj.iweight;
            md.nonlinprice = strcmp(obj.var.price, obj.nonlinparams);
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
                newmarket.p = obj.p(newmarket.selection, 1);
                newmarket.v = [];
                for k = 1:size(obj.x2, 2)
                    temp = obj.v(newmarket.selection,:,k);
                    newmarket.v = temp(1,:); 
                end
                newmarket.init();
                obj.period{t} = newmarket;
            end
        end
            
        function randdraws(obj, selection)
            K = size(obj.x2, 2); % Number of variables
            obj.v = [];
            if strcmpi(obj.settings.drawmethod, 'quadrature')
                [X, obj.iweight] = nwspgr('KPN', K, obj.settings.quaddraws);
                obj.settings.nind = length(obj.iweight);
                for k = 1:K
                    obj.v(:,:,k) = repmat(X(:,k)', size(obj.x2, 1), 1);
                end
            else
                obj.iweight = ones(obj.settings.nind, 1) / obj.settings.nind ;
                if isempty(obj.draws) % Can be set manually
                    if obj.settings.marketdraws
                        obj.draws = Sampling.draw(obj.settings.drawmethod, ...
                            K, max(obj.marketid) * obj.settings.nind, ...
                            obj.config.randstream);
                    else
                        obj.draws = Sampling.draw(obj.settings.drawmethod, ...
                            K, obj.settings.nind, obj.config.randstream);
                    end
                end
                for k = 1:K
                    % For each k, create draws for each market - wk (nind x nmkt)
                    if obj.settings.marketdraws
                        wk = reshape( obj.draws(:, k), obj.settings.nind, ...
                            max(obj.marketid))';
                    else
                        wk = repmat( obj.draws(:, k)', max(obj.marketid), 1);
                    end
                    % Duplicate draws for all products in each market
                    % Each k can have a different type of distribution
                    if obj.nonlintype(k) == 0
                        obj.v(:,:,k) = wk(obj.marketid, :);
                    elseif obj.nonlintype(k) == 1 % log-normal
                        obj.v(:,:,k) = exp(wk(obj.marketid, :)); 
                    elseif obj.nonlintype(k) == 2 % triangular
                        obj.v(:,:,k) = Sampling.triangular(wk(obj.marketid, :)); 
                    end              
                end
            end
            if ~isempty(selection)
                vnew = [];
                % This could be done above in creation, but this is simpler
                for k = 1:K
                    vnew(:,:,k) = obj.v(selection, :, k);
                end
                obj.v = vnew;
            end
        end
        
        function newedelta = finddelta(obj, rc_sigma, edelta, tolerance)
            sel = logical(obj.dummarket);
            newedelta = zeros(size(edelta));
            for t = 1:max(obj.marketid)
                newedelta(sel(:,t)) = obj.period{t}.finddelta(rc_sigma, ...
                    edelta(sel(:,t)), tolerance);
            end
        end       
    end
end

