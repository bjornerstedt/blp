classdef Merger < handle
    %MERGER Utility class to reduce coding of a merger
    %   Merger(demand) to create
    %   merger.merge(buyer, seller) to merge
    %   $Id: Merger.m 111 2015-04-30 11:18:10Z d3687-mb $
    
    properties
        demand
        market1
        market2
        selection
        results
        buyer
        seller
        var
    end
    
    methods        
        function result = merge(obj, demand)
            obj.demand = copy(demand);
            obj.demand.initSimulation(obj.selection);
            obj.market1 = Market(obj.demand);
            if ~isempty(obj.var.firm)
                obj.market1.var.firm = obj.var.firm;
            else
                error('Merger.firmvar has to be specified');
            end
            
            obj.market1.init();            
            obj.market1.findCosts();
            
            obj.market2 = copy(obj.market1);
            obj.market2.firm(obj.market2.firm == obj.seller ) = obj.buyer;
            obj.market2.p0 = obj.market1.p;
            obj.market2.init(); 
            flag = obj.market2.equilibrium(); 
            obj.results.isequilibrium = (flag > 0);
            
            [result, mpc] = obj.market1.compare(obj.market2.p);
            disp 'Merger results'
            disp(result)
            disp 'Average Price increase'
            disp(mpc) 
            obj.results.comparison = result;
        end
        
        function delete(obj)
            delete(obj.demand);
            delete(obj.market1);
            delete(obj.market2);
        end
    end
end

