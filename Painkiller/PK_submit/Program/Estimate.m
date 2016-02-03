classdef Estimate  < matlab.mixin.Copyable 
    % Basic table and logit functionality
    %   Subclassed by DemandEstimate, NestedLogit and MixedLogit
    %   $Id: Estimate.m 115 2015-05-21 09:11:46Z d3687-mb $
    
    properties
        settings % Structure with settings, defined in constructor
        config % Structure with internal settings, defined in constructor
        var % Structure with variable names
        panelid
        marketid
        nmkt
        Xorig % Non-demeaned vars
        Zorig % Non-demeaned vars
        y
        X
        Z
        results   
        beta  % Move to results?
        betadummies
    end
    properties % (SetAccess = protected )
        T
        vars  % vars to display in output      
        dummyvars = [];
        period % Cell array used for period / market calcs
        lsdv
    end
    
    methods

        function R = init(obj, varargin)
            % INIT: Initialize estimation of logit
            % 1) Generates nested logit log-shares
            % 2) Creates LSDV and constant or demeans all variables
            % 3) Adds dummy variables
            % 4) Creates X, Z and y matrices
            if nargin > 1
                selection = varargin{1};
                obj.T = obj.T(selection, :);
            else
                selection = [];
            end
            if ~isempty(obj.var.instruments)
                instrumentnames = strsplit(strtrim(obj.var.instruments));
            else
                instrumentnames = [];
            end
            exog = strsplit(strtrim(obj.var.exog));
            if ~isempty(obj.var.endog)
                endog = strsplit(strtrim(obj.var.endog));
            else
                endog = [];
            end
            if ~isempty(obj.dummyvars)
                exog = [exog,  obj.dummyvars.Properties.VariableNames];
            end
            if ~strcmpi(obj.settings.paneltype, 'none')
            [~,~,id] = unique(obj.T{:,strsplit(strtrim(obj.var.market))}, 'rows');
            obj.marketid = id;
            try
                obj.marketid(1,1); 
                % note that marketsize and quantity not tested as create
                % will create these.
            catch err
                display(err);
                error('Variables have to be defined')
            end
            obj.nmkt = max(obj.marketid);  
            end
            if ~isempty(obj.var.depvar)
                obj.y = obj.T{:, obj.var.depvar};
            end
            T = obj.T;
            varlist = [endog, exog, instrumentnames];
            if ~strcmpi(obj.settings.paneltype, 'none')
                T = T(:, [obj.var.panel, varlist] );
                panel = strsplit(strtrim(obj.var.panel));
                obj.panelid =  T{:, panel };
                obj.panelid(1,1);
            end
            constant = [];
            if strcmpi(obj.settings.paneltype, 'lsdv')
                if obj.settings.nocons 
                    obj.vars =[endog, exog];
                else
                    obj.vars =[endog, exog, { 'constant'}];
                    constant = ones(size(T, 1), 1);
                end
                if isempty(obj.lsdv) || isempty(selection)
                    obj.lsdv = dummyvar(grp2idx(T{:, panel } ));
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
            R = T;
            %Definitions of matrices
            obj.Xorig = [obj.T{: , exog} constant obj.lsdv ];
            obj.Zorig = [obj.T{: , [exog, instrumentnames]} constant obj.lsdv ];
            obj.X = [T{: , [endog, exog]} constant obj.lsdv ];
            obj.Z = [T{: , [exog, instrumentnames]} constant obj.lsdv  ];
            if any(isnan(table2array(T)))
                warning('NaN occurs in data')
            end
        end
        
        function f = addDummy(obj, name)
            dummy = array2table(dummyvar(obj.T{:,{name}}));
            dummy.Properties.VariableNames = ...
                arrayfun(@(x) {sprintf('%s%d', name, x)}, 1:size(dummy,2));
            dummy(:,1) = [];
            obj.dummyvars = [obj.dummyvars, dummy];
            f = dummy;
        end
        
        function createVarcovar(obj, varcovar)
            varsel = 1:length(obj.vars );
            varcovar = array2table(varcovar(varsel, varsel) );
            varcovar.Properties.VariableNames = obj.vars;
            varcovar.Properties.RowNames = obj.vars;
            obj.results.params.varcovar = varcovar;
        end
        
        function resultTables(obj)
            var = fieldnames(obj.var);
            value = struct2cell(obj.var);
            obj.results.var = table(var, value);
            setting = fieldnames(obj.settings);
            value = struct2cell(obj.settings);
            obj.results.settings = table(setting, value);    
            if ~isempty(obj.results.other)
                obj.results.properties = ...
                    table( struct2cell(obj.results.other) );
                obj.results.properties.Properties.RowNames = ...
                    fieldnames(obj.results.other);
            end
        end
            
        function R = createResults(obj, beta, varcovar, varnames)
            se = sqrt(  diag(varcovar) );  
            tvalue = beta ./ se;
            obj.results.params.beta = obj.beta;
            R = array2table([beta, se, tvalue ]);
            R.Properties.VariableNames = {'Coef','Std_err', 't_value'};
            R.Properties.RowNames = varnames;
            obj.results.estimate = R;
            obj.resultTables();
        end
        
        function R = estimate(obj)
            switch obj.settings.estimateMethod 
                case 'ols'
                    [beta, varcovar] = obj.ols();
                case 'gls' 
                    [beta, varcovar] = obj.gls();
                case 'gmm'
                    [beta, varcovar] = obj.gmm();
                otherwise
                    error(['estimateMethod: ', obj.settings.estimateMethod, ...
                        ' is not defined']);
            end
            obj.beta = beta;
            if strcmpi(obj.settings.paneltype, 'lsdv')
            varsel = 1:length(obj.vars );
                varcovar = varcovar(varsel, varsel) ;
                beta = beta(varsel, :);
                obj.results.betadummies = obj.beta(length(obj.vars)+1:end);
            end
            obj.results.other = [];
            R = obj.createResults(beta, varcovar, obj.vars);            
            obj.createVarcovar(varcovar);
        end
        
        function [beta, varcovar] = ols(obj)
            beta = inv(obj.X'*obj.X)*(obj.X'*obj.y);
            xi = obj.y - obj.X*beta;
            dgf = (size(obj.X,1) - size(obj.X,2));
            varcovar = ((xi'*xi) ./ dgf) * inv(obj.X'*obj.X);
        end
        
        function [beta, varcovar] = gls(obj)
            y = obj.y;
            X = obj.X;
            Z = obj.Z;
            dgf = (size(obj.X,1) - size(obj.X,2));
            
            W = inv(Z'*Z);
            mid = Z*W*Z';
            sst = inv(X'*mid*X);
            beta = sst * (X'*mid*y);
            xi = y - X * beta;
            % Add robust estimate
            ser = (xi'*xi) ./ dgf;
            varcovar = ser * sst;
        end
        
        function [beta, varcovar] = gmm(obj)
            y = obj.y;
            X = obj.X;
            Z = obj.Z;
            dgf = (size(obj.X,1) - size(obj.X,2));
            
            %STAGE I: INITIAL WEIGHTING MATRIX*/
            W = inv(Z'*Z);
            mid = Z*W*Z';
            btsls = inv(X'*mid*X) * (X'*mid*y);
            xi = y - X*btsls;

            %STAGE II: OPTIMAL WEIGHTING MATRIX
            W = inv((bsxfun(@times,xi,Z))'*(bsxfun(@times,xi,Z)));
            mid = Z*W*Z';
            beta = inv(X'*mid*X)*(X'*mid*y);
            xi = y - X*beta;
            varcovar = inv(X'*mid*X);
        end

        % Set settings or config parameters with struct or cell array of
        % struct
        function fexist = set(obj, prefstruct)
            fexist = true;
            if iscell(prefstruct)
                for i = 1:length(prefstruct)
                    fexist = fexist && obj.set(prefstruct{i});
                end
            else
                names = fieldnames(prefstruct)
                for i = 1:length(names)
                    if isfield(obj.settings, names{i}) 
                        obj.settings.(names{i}) = prefstruct.(names{i});
                    elseif isfield(obj.config, names{i}) 
                        obj.config.(names{i}) = prefstruct.(names{i});
                    else
                        fexist = false;
                    end            
                end 
            end
        end
        
        function obj = Estimate(varargin) 
            if nargin > 0
                obj.T = varargin{1};
            end
            obj.config.quietly = false;
            obj.config.isdemeaned = false; 
            obj.config.estimateDescription = 'Linear Estimate'; 
            obj.settings.paneltype = 'fe';
            obj.settings.nocons = false;   
            obj.settings.estimateMethod = 'ols';
            
            obj.var.market = ''; 
            obj.var.panel = '';
            obj.var.depvar = '';
            % Sets of variables
            obj.var.exog = ''; 
            obj.var.endog = '';
            obj.var.instruments = '';
        end        
        
        function f  = demean(obj, T, vars, indexvars )
        %DEMEAN Subtract average value of vars
        %   The index variables are the variables averages are taken over.
            [uniqueRows,~,rowIdx]=unique(obj.T{:,indexvars},'rows');
            if length(vars) > 0
                T = T( :, [indexvars vars]);
            else
                vars = T.Properties.VariableNames;
            end
            B = varfun(@mean, T,'GroupingVariables', indexvars);
            B =  B{rowIdx, length(indexvars)+2:end};
            T(: , indexvars) = []; 
            B = array2table(table2array(T) - B );
            B.Properties.VariableNames = vars;
            f = B;
        end
    end
end

