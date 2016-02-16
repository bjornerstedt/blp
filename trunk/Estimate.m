classdef Estimate  < matlab.mixin.Copyable 
    % Basic table and logit functionality
    %   Subclassed by NestedLogit and MixedLogit
    
    properties
        data % Data table used by init
        panelid
        results   
        y
        X
        Z
        beta  % Move to results?
        settings % Structure with settings, defined in constructor
        config % Structure with internal settings, defined in constructor
        var % Structure with variable names
        Xorig % Non-demeaned vars
        Zorig % Non-demeaned vars
    end
    properties (SetAccess = protected, Hidden = true )
        vars  % vars to display in output      
        period % Cell array used for period / market calcs
        lsdv
    end
    
    methods
        function R = estimate(obj, varargin)
            obj.init(varargin{:});
            switch lower(obj.settings.estimateMethod)
                case 'ols'
                    [beta, varcovar] = obj.ols();
                case '2sls' % Change name to 2sls Make default if obj.Z~=[]
                    [beta, varcovar] = obj.tsls();
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
                obj.results.params.betadummies = obj.beta(length(obj.vars)+1:end);
            end
%             obj.results.other = [];
            R = obj.createResults(beta, varcovar, obj.vars);            
            obj.createVarcovar(varcovar);
        end
          
        function fexist = set(obj, prefstruct)
        % Set settings or config parameters with struct or cell array of
        % struct
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
                obj.data = varargin{1};
            end
            obj.var = SettingsClass({'market','panel','depvar','exog', ...
                'endog','instruments'});
            obj.settings = SettingsClass({'paneltype', 'nocons', ...
                'estimateMethod', 'robust', 'logLinear'});
%            obj.config = SettingsClass({'quietly','nocons ss','estimateMethod'});
            
            obj.config.quietly = true;
            obj.results.estimateDescription = 'Linear Estimate'; 
            
            obj.settings.paneltype = 'fe';
            obj.settings.nocons = false;   
            obj.settings.estimateMethod = 'ols';
            obj.settings.robust = true;
        end
        
       function selection = init(obj, varargin)
            % init([selection]): Initialize estimation of logit, with selection as
            % optional argument. 
            % Creates X, y and Z, as well as panelid
            args = inputParser;
            args.addParameter('selection', [],  @islogical );
            args.parse(varargin{:});
            selection = args.Results.selection;
            exog = [];
            endog = [];
            instruments = [];
            if ~isempty(obj.var.instruments)
                instruments = strsplit(strtrim(obj.var.instruments));
            else
                obj.settings.estimateMethod = 'ols';
            end
            if ~isempty(obj.var.exog)
                exog = strsplit(strtrim(obj.var.exog));
            end
            if ~isempty(obj.var.endog)
                endog = strsplit(strtrim(obj.var.endog));
            end
            % initAdditional creates Xorig and returns names of created
            [obj.Xorig, created] = obj.initAdditional( );
            if ~strcmpi(obj.settings.paneltype, 'none')
                if isempty(obj.var.panel) || ~obj.isvar(obj.var.panel, obj.data)
                    error('var.panel has to be specified')
                end
                panel = strsplit(strtrim(obj.var.panel));
                [~, ~,obj.panelid] =  unique(obj.data{:, panel }); 
            end
            % Can have log model here.
            if obj.settings.nocons || strcmpi(obj.settings.paneltype, 'fe')
                obj.vars = [created, endog, exog];
                constant = [];
            else
                obj.vars =[created, endog, exog, { 'constant'}];
                constant = ones(size(obj.data, 1), 1);
            end
            switch lower(obj.settings.paneltype)
                case 'lsdv'
                    if isempty(obj.lsdv) || isempty(selection)
                        obj.lsdv = dummyvar(grp2idx(obj.data{:, panel } ));
                        obj.lsdv = obj.lsdv(:, sum(obj.lsdv)~=0 );
                        if ~obj.settings.nocons
                            obj.lsdv = obj.lsdv(:, 2:end);
                        end
                    else
                        % obj.lsdv cannot be recreated on a selection as it
                        % screws up dummy count and numbering
                    end
                case 'fe'
                    obj.lsdv = [];
                case 'none'
                    obj.lsdv = [];
                otherwise
                    error('Unknown paneltype')
            end
            if ~isempty(obj.var.depvar)
                obj.y = obj.data{:, obj.var.depvar};
            end
            obj.Xorig = [obj.Xorig, obj.data{:, [endog, exog]}, constant, obj.lsdv];
            if ~isempty(obj.var.instruments)
                obj.Zorig = [obj.data{:, [exog, instruments]}, constant, obj.lsdv];
            end
            if strcmpi(obj.settings.paneltype, 'fe')
                obj.X = obj.demean(obj.Xorig);
                obj.y = obj.demean(obj.y);
                obj.Z =  obj.demean(obj.Zorig);
                obj.results.dgf = (size(obj.X,1) - size(obj.X,2)) ...
                    - length(unique(obj.panelid));
            else
                obj.X = obj.Xorig;
                obj.Z =  obj.Zorig;
                obj.results.dgf = (size(obj.X,1) - size(obj.X,2));
            end
            obj.results.observations = size(obj.X,1);
            if ~isempty(selection)                
                obj.X = obj.X(selection, :);
                obj.y = obj.y(selection);
                obj.Z = obj.Z(selection, :);
                obj.Xorig = obj.Xorig(selection, :);
                obj.Zorig = obj.Zorig(selection, :);
            end
%             if any(isnan(table2array(T)))
%                 warning('NaN occurs in data')
%             end
        end
        
    end
    
    methods(Access = protected, Hidden = true)
        
        function [Xorig, created] = initAdditional(obj)
            % Additional initialisation in subclasses
            created = [];
            Xorig = [];
        end
               
        function createVarcovar(obj, varcovar)
            varsel = 1:length(obj.vars );
            varcovar = array2table(varcovar(varsel, varsel) );
            varcovar.Properties.VariableNames = obj.vars;
            varcovar.Properties.RowNames = obj.vars;
            obj.results.params.varcovar = varcovar;
        end
        
        function resultTables(obj)
%             var = fieldnames(obj.var);
%             value = struct2cell(obj.var);
%             obj.results.var = table(var, value);
            obj.results.var = obj.var.getPropertyTable();
            obj.results.settings = obj.settings.getPropertyTable();  
%             if ~isempty(obj.results.other)
%                 obj.results.properties = ...
%                     table( struct2cell(obj.results.other) );
%                 obj.results.properties.Properties.RowNames = ...
%                     fieldnames(obj.results.other);
%             end
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
        
       % Put beta and varcovar in obj only in estimate()
        function [beta, varcovar] = ols(obj)
            invXX = inv(obj.X' * obj.X);
            beta = invXX * (obj.X' * obj.y);
            xi = obj.y - obj.X * beta;
            if obj.settings.robust
                xiX =bsxfun(@times,xi,obj.X);
                varcovar = invXX * (xiX' * xiX) * invXX;
            else
                varcovar = ((xi' * xi) ./ obj.results.dgf) * inv(obj.X' * obj.X);
            end
        end
        
        function [beta, varcovar] = tsls(obj)
            W = inv(obj.Z' * obj.Z);
            Zprod = obj.Z * W * obj.Z';
            invZZ = inv(obj.X' * Zprod * obj.X);
            beta = invZZ * (obj.X' * Zprod * obj.y);
            
            % See docs/design.lyx, section Robust errors
            xi = obj.y - obj.X * beta;
            if obj.settings.robust
                xiXhat = bsxfun(@times, xi, Zprod * obj.X);
                varcovar = invZZ * (xiXhat' * xiXhat) * invZZ;
                % TODO: Add dgf adjustment
            else
                varcovar = ((xi' * xi) ./ obj.results.dgf) * invZZ;
            end
        end
        
        function [beta, varcovar] = gmm(obj)            
            % STAGE I: INITIAL WEIGHTING MATRIX
            W = inv(obj.Z' * obj.Z);
            mid = obj.Z * W * obj.Z';
            btsls = inv(obj.X' * mid * obj.X) * (obj.X' * mid * obj.y);
            xi = obj.y - obj.X * btsls;

            % STAGE II: OPTIMAL WEIGHTING MATRIX
            xiZ = bsxfun(@times,xi,obj.Z);
            W = inv(xiZ' * xiZ);
            mid = obj.Z * W * obj.Z';
            beta = inv(obj.X' * mid * obj.X) * (obj.X' * mid * obj.y);
            xi = obj.y - obj.X * beta;

            % See docs/design.lyx, section Robust errors
            % Note that second stage W is not used.
            varcovar = inv(obj.X' * mid * obj.X);
            % TODO: Add dgf adjustment
        end
       
        function mn  = demean(obj, data )
        %DEMEAN Subtract the mean value of vars over panel
            if istable(data)
                data = table2array(data);
            end
            if isempty(data)
                mn = [];
                return
            end
            mn = zeros(max(obj.panelid), size(data, 2));
            for i = 1:size(data, 2)
                mn(:, i) = accumarray(obj.panelid, data(:, i), [], @mean );
            end
            mn =  mn(obj.panelid, :);
            mn = data - mn;
        end
        
        function cpObj = copyElement(obj)
        % Overrides copyElement method
        
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Make a deep copy of the DeepCp object
            cpObj.var = copy(obj.var);
            cpObj.settings = copy(obj.settings);
            if isa(obj.config, 'SettingsClass')
                cpObj.config = copy(obj.config);
            end
        end
                
    end
    
    methods( Static )
        function [instruments, names] = countInstruments(data, market, typeNames, varargin)
        % countInstruments(market, typeNames, sumVars, dataTable) 
        % sums vars by market over firm and typenames 
        % Market has to be specified, as the function is invoked prior to
        % the creation of marketid in init.
            instruments = table;
            k = 1;
            if nargin == 4
                sumVar = data{ :, varargin{1}};                
            else
                sumVar = ones(size(data,1), 1);
            end
            names = cell(0);
            for s = 1:size(sumVar, 2);
                for i = 0:length(typeNames)
                    sels = nchoosek(1:length(typeNames), i);
                    for j = 1:size(sels, 1)
                        sel = [market, typeNames(sels(j,:))];
                        names = [names; {sel}];
                        [~, ~, index] = unique(data(:, sel));
                        count = accumarray(index, sumVar(:, s) );
                        instruments.(sprintf('inst%d',k)) = count(index) - sumVar(:, s);
                        k = k + 1;
                    end
                end
            end
            instruments = table2array(instruments);
        end
    end
    
    methods(Static , Hidden = true )
        function v = isvar(x,y)
            v = any(strcmp(x, y.Properties.VariableNames));
        end
        
        function R = mean(index, data)
            [~,~, index] = unique(index);
            R = accumarray(index, data);
        end
     
        function R = means(data, cols, index, varargin)
        % means(data, tableVars, indexVar, weightVar) calculates means of cols 
        % over and index. With optional weight variable, weighted means are 
        % calculated. Weights can either be specified as the variable name 
        % in the table data or as a vector.
            [cats, ~, rowIdx] = unique(data( : , index), 'rows');
            if ~isempty(varargin) && ~isempty(varargin{1}) 
                if ischar(varargin{1})
                    wt0 = data.(varargin{1});
                else
                    wt0 = varargin{1};
                end
                wt = accumarray(rowIdx, wt0);
                wt = wt0 ./ wt(rowIdx,:);
                R = array2table(splitapply(@(x)sum(x,1), bsxfun(@times, data{:, cols}, ...
                    wt), rowIdx), 'VariableNames', cols);
            else
                R = array2table(splitapply(@(x)mean(x,1), data{:, cols}, ...
                    rowIdx), 'VariableNames', cols);
            end
            R = [cats, R];
        end
    end
end

