classdef ExcelTable < handle
    %ExcelTable Write tables to Excel
    %   Keeps track of positions of tables in sheet
    %   $Id: ExcelTable.m 112 2015-05-07 16:05:37Z d3687-mb $
    
    properties
        rowmargin = 2 % Empty rows between tables + heading size
        colmargin = 2
    end
    properties (SetAccess = private )
        filename
        sheet
        row = 1
        col = 1
        height = 0 
        width = 0 
        basecol = 1
        colwise = 0
    end
    
    methods
        function setSheet(obj, sheet)
            obj.sheet = sheet;
            obj.row = 1;
            obj.col = 1;
            obj.width = 0;
            obj.height = 0;
            obj.basecol = 1;
            obj.colwise = 0;
        end
        
        function setColwise(obj, col)
            if isstr(col)
                if length(col) > 1
                    error('setCol: col must be number or single char A-Z');
                end
                col = double(col) - 64;
            end
            obj.basecol = col;
            obj.col = col;
            if obj.colwise == 0
                obj.colwise = obj.row;
            else
                obj.row = obj.colwise;
            end
            obj.width = 0;
            obj.height = 0;
        end
        
        function  write(obj, outtable, varargin)
            p = inputParser;
            p.addRequired('outtable', @istable);
            p.addParameter('sheet', []);
            p.addParameter('heading', [], @isstr);
            p.addParameter('rownamesheading', [], @isstr);
            p.addParameter('below', true, @islogical);
            p.parse(outtable, varargin{:});
            if ~isempty(p.Results.sheet)
                obj.sheet =  p.Results.sheet;
            end
            options = {'WriteRowNames', true};
            if ~isempty(obj.sheet)
                options = [options, {'sheet', obj.sheet}];
            end
            obj.width = size(outtable,2) + obj.rowmargin;
            if p.Results.below
                obj.row = obj.row + obj.height;
                obj.col = obj.basecol;
                obj.height = size(outtable,1) + obj.rowmargin;
            else
                obj.col = obj.col + obj.width;
                obj.height = max(obj.height, size(outtable,1)+obj.rowmargin);
            end
            if ~isempty(p.Results.heading)
                xlswrite(obj.filename, {p.Results.heading}, obj.sheet, ...
                    ExcelTable.makerc(obj.row, obj.col) );
                hm = 2;
                obj.height = obj.height + hm;
            else
                hm = 0;
            end
            rc = ExcelTable.makerc(obj.row + hm, obj.col);
            writetable(outtable, obj.filename, 'range', rc, options{:});
            
            if ~isempty(p.Results.rownamesheading)
                xlswrite(obj.filename, {p.Results.rownamesheading}, ...
                    obj.sheet, ExcelTable.makerc(obj.row + hm, obj.col) );
            end

        end
        
        function obj = ExcelTable(name, varargin)
            p = inputParser;
            p.addRequired('name', @isstr);
            p.addParameter('sheet',[], @isstr);
            if nargin > 1
            p.parse(name, varargin{:});
            obj.sheet =  p.Results.sheet;
            end
            obj.filename = name;
        end
        
    end
    
    methods(Static)
        function rc = makerc(row, col)
            if col <= 26
                z = char(mod(col-1, 26)+65);
            else
                z = [char(floor((col-1)/26)+64), char(mod(col-1, 26)+65)];
            end
            rc = sprintf('%s%d', z, row);
        end
    end
end

