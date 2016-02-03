classdef Demand
    %DEMAND Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n = 100
    end
    properties (SetAccess = protected)
        
    end
    
    methods
        function f = parameters(obj, x)
            f = obj.n^2+x;
        end
        
        function f = objective(obj, x)
            f = (x-obj.n)^2
        end
        function f = estimate(obj, x)
            %options = optimset('Display','iter-detailed','TolX',1e-6,'TolFun',1e-6);
            %[rc_sigma,fval,exitflag,output] = fminsearch(@(x)objective(x,data),rc_sigma0);
            [y,fval,exitflag,output] = fminsearch(@obj.objective, 6);
            f= y
        end
    end
end

