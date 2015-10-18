classdef DemandEstimate  < Estimate 
    % Basic table and logit functionality
    %   Subclassed by NestedLogit and MixedLogit
    %   $Id: DemandEstimate.m 118 2015-05-26 10:34:48Z d3687-mb $
    
    properties        
        nestlist
        s
        s0
        ls
        p
        q 
        ms
        dummarket 
    end
        
    methods
        function name = getPriceName(obj)
            if obj.settings.ces
                name = 'lP';
            else
                name = obj.var.price;
            end
        end
        
        function names = getLogShareNames(obj)
            lsnames = {[],{'lsjg'},{'lsjh', 'lshg'}};
            names = lsnames{length(obj.nestlist)+1};
        end
        
        function price = getPrice(obj)
            if obj.settings.ces
                price = log(obj.T{:,obj.var.price});
            else
                price = obj.T{:,obj.var.price};
            end
        end

        function init(obj, varargin)
            if nargin > 1
                selection = varargin{1};
                obj.T = obj.T(selection, :);
            else
                selection = [];
            end
            if ~isempty(obj.var.nests)
                obj.nestlist = strsplit(strtrim(obj.var.nests));
            end
            endog = [obj.getPriceName(), obj.getLogShareNames(), obj.var.endog];
            exog = strsplit(strtrim(obj.var.exog));
            if ~isempty(obj.dummyvars)
                exog = [exog,  obj.dummyvars.Properties.VariableNames];
            end
            [~,~,id] = unique(obj.T{:,strsplit(strtrim(obj.var.market))}, 'rows');
            obj.marketid = id;
            obj.dummarket = dummyvar(obj.marketid);
            obj.nmkt = max(obj.marketid);   % number of markets
            obj.p = obj.T{:, obj.var.price}; % Used in simulation 
            obj.ms = obj.T{: , obj.var.marketsize};
            try
                obj.p(1,1); obj.marketid(1,1); 
            catch err
                display(err);
                error('Variables have to be defined')
            end
            if ~isempty(obj.var.instruments)
                instrumentnames = strsplit(strtrim(obj.var.instruments));
            else
                instrumentnames = [];
            end
             if ~isempty(obj.var.quantity) 
                obj.q = obj.T{: , obj.var.quantity};
                sharenames ={'s', 's0'};
                shares = obj.generateShares(); 
                varlist = [obj.var.depvar, endog, exog, instrumentnames, sharenames];
                T = [obj.T shares obj.dummyvars];
            else  
                T = obj.T;
                varlist = [endog, exog, instrumentnames];
            end
            if obj.settings.ces
                T.lP = log(obj.T{:,obj.var.price});
            end
            if ~strcmpi(obj.settings.paneltype, 'none')
                T = T(:, [obj.var.panel, varlist] );
                panel = strsplit(strtrim(obj.var.panel));
                obj.panelid =  T{:, panel };
                obj.panelid(1,1);
            end
            if ~isempty(obj.var.quantity)
                obj.s = T.s;
                obj.s0 = T.s0;
                obj.ls = T.ls;
            end
            constant = [];
            if strcmpi(obj.settings.paneltype, 'lsdv')
                if obj.settings.nocons 
                    obj.vars =[endog, exog];
                else
                    obj.vars =[endog, exog, { 'constant'}];
                    constant = ones(size(obj.T, 1), 1);
                end
                if isempty(obj.lsdv) || isempty(selection)
                    obj.lsdv = dummyvar(grp2idx(obj.T{:, panel } ));
                    obj.lsdv = obj.lsdv(:, sum(obj.lsdv)~=0 );
                    if ~obj.settings.nocons
                        obj.lsdv = obj.lsdv(:, 2:end);
                    end
                else
                    obj.lsdv = obj.lsdv(selection,:);
                end
            elseif strcmpi(obj.settings.paneltype, 'fe')
                obj.vars = [endog, exog];
                if ~obj.config.isdemeaned
                    T = obj.demean(T, varlist, panel);
                end
            elseif strcmpi(obj.settings.paneltype, 'none')
                if obj.settings.nocons 
                    obj.vars =[endog, exog];
                else
                    obj.vars =[endog, exog, { 'constant'}];
                    constant = ones(size(T, 1), 1);
                end
            else
                error('Unknown paneltype')
            end
            obj.X = [T{: , [endog, exog]} constant obj.lsdv ];
            if ~isempty(instrumentnames)  
                obj.Z = [T{: , [exog, instrumentnames]} constant obj.lsdv  ];
            end
            if ~isempty(obj.var.quantity) && ~isempty(obj.var.depvar) % HACK!!!
                obj.y = T{:, obj.var.depvar};
            end
            if any(isnan(table2array(T)))
                warning('NaN occurs in data')
            end
        end
        
        function resultTables(obj)
            resultTables@Estimate(obj);            
            if ~obj.config.quietly
                disp(['Estimate of: ', obj.config.estimateDescription])
                disp(obj.results.estimate);
                s0 = 1 - obj.dummarket' * obj.s;
                disp('Share of outside good');
                fprintf('  Mean: %0.3f Min: %0.3f Max: %0.3f \n', ...
                    mean(s0), min(s0), max(s0));
            end            
        end
        
        function obj = DemandEstimate(varargin)
            if nargin > 0
                obj.T = varargin{1};
            end
            obj.var.depvar = 'ls';
            obj.var.quantity = '';
            obj.var.price = '';        
            obj.var.nests = '';
            obj.var.marketsize = '';
            obj.settings.ces = false;
        end        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods % (Access = protected)
        function S = generateShares(obj )
            %SUBTOTALS Sums sumvar in table T by index category variable list
            %   Mimics Stata by index: egen subtotal = total(sumvar)
            if obj.settings.ces
                Q = obj.q .* obj.p;
            else
                Q = obj.q;
            end
            subtotal = accumarray(obj.marketid, Q);
            subtotal =  subtotal(obj.marketid,:);
            S = table();
            S.s = Q ./ obj.ms;
            S.s0 = 1 - ( subtotal ./ obj.ms );
            S.ls = log(S.s ./ S.s0);
            market = strsplit(strtrim(obj.var.market));
            if length(obj.nestlist) >= 1
                groupsubtotal =  obj.subtotals(Q, ...
                    [market obj.nestlist(1)]);
                S.lsjg = log( Q ./ groupsubtotal);
            end
            if length(obj.nestlist) == 2
                subgroupsubtotal =  obj.subtotals(Q, ...
                    [market obj.nestlist]);
                S.lsjh = log( Q ./ subgroupsubtotal);
                S.lshg = log( subgroupsubtotal ./ groupsubtotal );
            end
        end
        
        function f = subtotals(obj, sumvar, index )
            %SUBTOTALS Sums sumvar in table T by index category variable list
            %   Mimics Stata by index: egen subtotal = total(sumvar)
                [~, ~, rowIdx] = unique(obj.T( : , index), 'rows');
                subtotal = accumarray(rowIdx, sumvar);
                f =  subtotal(rowIdx,:);
        end
        
        % Standard error of the regression
        % Cameron & Trivedi p 287
        function sd = sdreg(obj)
            sm = obj.s - mean(obj.s);
            sp = obj.shares(obj.p, 1); % Compute predicted shares
            e = obj.s - sp;
            sd = sqrt(e'*e ./(size(obj.X,1)-length(obj.vars2)-length(obj.vars2)));
        end        
    end
end

