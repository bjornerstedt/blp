classdef SettingsClass < dynamicprops & matlab.mixin.Copyable 
    %Settings Handle settings in a better way than with a struct 
    %   Initialize with SettingsClass(x) where x is a string or cell array
    %   Set values from struct s with init(x)
    %   Get struct with properties wit getProperties()
    %   Additional properties can be set with set(x)
    %   A common problem with using a struct of settings that can be set by
    %   the user is that the set of properties is not fixed. Misspelling a
    %   property leads to a new field being created.
   
    methods
        function setParameters(obj, setnames)
            if ischar(setnames)
                obj.addprop(setnames);
            else
                if ~iscell(setnames)
                    error('SettingsClass can only be set with string or cell array');
                end
                for i=1:length(setnames)
                    obj.addprop(setnames{i});
                end
            end
        end
        
        function init(obj, valstruct)
            if ~isstruct(valstruct)
                    error('init can only be invoked with a struct');
            end
            names = fieldnames(valstruct);
            for i = 1:length(names)
                obj.(names{i}) = valstruct.(names{i});
            end
        end
        
        function str = getProperties(obj)
            names = fieldnames(obj);
            str = struct();
            for i = 1:length(names)
                str.(names{i}) = obj.(names{i});
            end
        end
        
        function tab = getPropertyTable(obj)
            names = fieldnames(obj);
            values = cell(size(names));
            for i = 1:length(names)
                values{i} = obj.(names{i});
            end
            tab = table(names, values);
        end
                
        function obj = SettingsClass(setnames)
            obj.setParameters(setnames);
        end
    end
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Make a deep copy of the DeepCp object
            names = fieldnames(obj);
            for i = 1:length(names)
                cpObj.addprop(names{i});
                cpObj.(names{i}) = obj.(names{i});
            end
        end
    end
end

