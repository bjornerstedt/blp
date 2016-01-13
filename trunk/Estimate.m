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
        dummyvars = []
        period % Cell array used for period / market calcs
        lsdv
    end
    
    methods
        function R = estimate(obj, varargin)
            obj.init(varargin{:});
            switch obj.settings.estimateMethod
                case 'ols'
                    [beta, varcovar] = obj.ols();
                case 'gls' % Change name to 2sls Make default if obj.Z~=[]
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
                'estimateMethod', 'robust'});
%            obj.config = SettingsClass({'quietly','nocons','estimateMethod'});
            
            obj.config.quietly = true;
            obj.results.estimateDescription = 'Linear Estimate'; 
            
            obj.settings.paneltype = 'fe';
            obj.settings.nocons = false;   
            obj.settings.estimateMethod = 'ols';
            obj.settings.robust = true;
        end
        
       function selection = init(obj, varargin)
            % INIT: Initialize estimation of logit, with selection as
            % optional argument. Create a always true selection to reset.
            args = inputParser;
            args.addParameter('selection', [],  @islogical );
            args.parse(varargin{:});
            selection = args.Results.selection;
            if ~isempty(obj.var.instruments)
                names.instruments = strsplit(strtrim(obj.var.instruments));
            else
                names.instruments = [];
                obj.settings.estimateMethod = 'ols';
            end
            names.exog = strsplit(strtrim(obj.var.exog));
            if ~isempty(obj.var.endog)
                names.endog = strsplit(strtrim(obj.var.endog));
            else
                names.endog = [];
            end
            if ~isempty(obj.dummyvars)
                names.exog = [names.exog,  obj.dummyvars.Properties.VariableNames];
            end
            
            names = obj.initAdditional( names, selection);            
            varlist = [obj.var.depvar, names.endog, names.exog, names.instruments];
           
            % Handle panel
            if ~strcmpi(obj.settings.paneltype, 'none')
                if isempty(obj.var.panel) || ~obj.isvar(obj.var.panel, obj.data)
                    error('var.panel has to be specified')
                end
                T = obj.data(:, [obj.var.panel, varlist] );
                panel = strsplit(strtrim(obj.var.panel));
                [~, ~,obj.panelid] =  unique(T{:, panel }); % Panelid should be created better
            else
                T = obj.data(:, varlist);             
            end
            if obj.settings.nocons || strcmpi(obj.settings.paneltype, 'fe')
                obj.vars =[names.endog, names.exog];
                constant = [];
            else
                obj.vars =[names.endog, names.exog, { 'constant'}];
                constant = ones(size(T, 1), 1);
            end
            switch lower(obj.settings.paneltype)
                case 'lsdv'
                    if isempty(obj.lsdv) || isempty(selection)
                        obj.lsdv = dummyvar(grp2idx(T{:, panel } ));
                        obj.lsdv = obj.lsdv(:, sum(obj.lsdv)~=0 );
                        if ~obj.settings.nocons
                            obj.lsdv = obj.lsdv(:, 2:end);
                        end
                    else
                        % obj.lsdv cannot be recreated on a selection as it
                        % screws up dummy count and numbering
                        %                obj.lsdv = obj.lsdv(selection,:);
                    end
                case 'fe'
                    T = obj.demean(T, varlist, panel);
                    obj.lsdv = [];
                case 'none'
                    obj.lsdv = [];
                otherwise
                    error('Unknown paneltype')
            end
            if ~isempty(obj.var.depvar)
                obj.y = T{:, obj.var.depvar};
            end
            newvars = [constant, obj.lsdv, obj.dummyvars];
            if ~isempty(selection)
                newvars = newvars(selection, :);
                T = T(selection, :);
                Torig = obj.data(selection, :);
            else
                Torig = obj.data;
            end
            obj.Xorig = [Torig{: , [names.endog, names.exog]}, newvars];
            obj.X = [T{: , [names.endog, names.exog]}, newvars];
            if ~isempty(obj.var.instruments)
                obj.Zorig = [Torig{: , [names.exog, names.instruments]}, newvars];
                obj.Z = [T{: , [names.exog, names.instruments]}, newvars];
            end
            if any(isnan(table2array(T)))
                warning('NaN occurs in data')
            end
        end
        
    end
    
    methods(Access = protected, Hidden = true)
        
        function names = initAdditional(obj, names, selection)
            % Additional initialisation in subclasses
        end
               
        function f = addDummy(obj, name)
            dummy = array2table(dummyvar(obj.data{:,{name}}));
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
            invXX = inv(obj.X'*obj.X);
            beta = invXX * (obj.X'*obj.y);
            xi = obj.y - obj.X*beta;
            if obj.settings.robust
                xiX =bsxfun(@times,xi,obj.X);
                varcovar = invXX*(xiX'*xiX)*invXX;
            else
                if strcmpi(obj.settings.paneltype, 'fe')
                    dgf = (size(obj.X,1) - size(obj.X,2)) - max(obj.panelid);
                else
                    dgf = (size(obj.X,1) - size(obj.X,2));
                end
                varcovar = ((xi'*xi) ./ dgf) * inv(obj.X'*obj.X);
            end
        end
        
        function [beta, varcovar] = gls(obj)
            y = obj.y;
            X = obj.X;
            Z = obj.Z;
            if strcmpi( obj.settings.paneltype, 'fe') 
                dgf = (size(obj.X,1) - size(obj.X,2)) - max(obj.panelid);
            else
                dgf = (size(obj.X,1) - size(obj.X,2));
            end
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
            
            % STAGE I: INITIAL WEIGHTING MATRIX
            W = inv(Z'*Z);
            mid = Z*W*Z';
            btsls = inv(X'*mid*X) * (X'*mid*y);
            xi = y - X*btsls;

            % STAGE II: OPTIMAL WEIGHTING MATRIX
            W = inv((bsxfun(@times,xi,Z))'*(bsxfun(@times,xi,Z)));
            mid = Z*W*Z';
            beta = inv(X'*mid*X)*(X'*mid*y);
            xi = y - X*beta;
            varcovar = inv(X'*mid*X);
        end
       
        function f  = demean(obj, T, vars, indexvars )
        %DEMEAN Subtract average value of vars
        %   The index variables are the variables averages are taken over.
            [uniqueRows,~,rowIdx]=unique(obj.data{:,indexvars},'rows');
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
    
    methods(Static , Hidden = true )
        function v = isvar(x,y)
            v = any(strcmp(x, y.Properties.VariableNames));
        end
        
        function R = means(data, cols, index, varargin)
        % means calculates means of cols over index. With optional weight
        % variable, weighted means are calculated. Weights can either be
        % specified as the variable name in the table data or as a vector.
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

