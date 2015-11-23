      % panel calcs for xi in MixedLogitDemand
      
function R = computeResults(obj, delta)
            [beta, xi] = obj.lpart(delta);
            obj.beta = beta;
            obj.alpha = obj.beta(strcmp(obj.getPriceName(), obj.vars));

            if strcmpi(obj.settings.paneltype, 'fe')
                index = obj.panelid;
                xa = [];
                for i = 1:size(obj.X,2)
                    xa = [xa accumarray(index, obj.X(:,i),[],@mean)];
                end
                xavt = xa(index, :);
                davt = accumarray(index, delta, [], @mean);
                davt = davt(index, :);
                muest = davt - xavt*beta;
                xi = delta - obj.X*beta - muest;
                obj.xi = xi;
            else
                obj.xi = xi;
            end
            
            %etc...
end
        % redefined optimalInstrument based on demeaned panel data
% 5. Compute the jacobian of delta wrt sigma
                obj.optimalInstrument = obj.deltaJacobian(obj.rc_sigma, exp(deltaHat));
                
% 6. Create optimal instruments
                % Panel correction, demean deltaJac
%                if obj.ispanel && ~obj.isLsdv 
                if strcmpi(obj.settings.paneltype, 'fe')
                    da = [];
                    for i = 1:size(obj.optimalInstrument ,2)
                        da = [da, accumarray(obj.panelid, obj.optimalInstrument(:, i),...
                            [],@mean)];
                    end
                    davt = da(obj.panelid, :);
                    obj.optimalInstrument  = (obj.optimalInstrument  - davt);
                end