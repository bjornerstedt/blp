classdef GmmMarket < Market
    % GmmMarket Adds simultaneous demand and cost estimation to Market class
    
    properties
        W
        ZWZ
        alphahist = []
    end
    
    methods
        function [beta, betaM] = minimize(obj)
            XWX = obj.X* inv(obj.X' * obj.X)* obj.X';
            W = inv(obj.D.Z' * obj.D.Z);
            obj.ZWZ = obj.D.Z * W * obj.D.Z';
            beta = obj.D.beta;
            dlen = length(obj.D.beta);
            if true
                options = optimset('Display','iter-detailed','MaxIter',50);
                options = optimset('Display','iter-detailed');
                [beta,fval,exitflag,output] = fminsearch( ...
                    @objective, beta, options);
            else
                options = optimoptions(@fminunc, 'Display','iter-detailed');
                [beta,fval,exitflag, output] = fminunc( ...
                    @objective, beta, options);
            end
            betaM = obj.beta;
          
            function dist = objective(beta)
                obj.D.setBeta(beta);
                obj.alphahist = [obj.alphahist ; obj.D.alpha];
                obj.findCosts();
                if min(obj.c) > 0
                    lc = log( obj.c );
                else
                    warning('Calculated negative costs');
                    lc = log( max(obj.c, .0001) );
                end
                obj.y = lc;
                obj.ols(); % Generate new beta
                xiD = obj.D.y - obj.D.X*beta;
                xiM = obj.y - obj.X*obj.beta;
                dist1 = xiD'* obj.ZWZ *xiD ;
                dist2 = xiM' *xiM;
                fprintf('Distance: %f, %f\n', dist1, dist2);
                dist = xiD'* obj.ZWZ *xiD + xiM'* XWX *xiM;
            end
            
        end
        
        function dist = objective2(obj, beta)
            e = obj.y - obj.X*beta;
            dist = e' * e;
        end
        
        function obj = GmmMarket(data)
            obj = obj@Market(data);
        end
    end
end

