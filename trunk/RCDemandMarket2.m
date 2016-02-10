classdef RCDemandMarket2  < matlab.mixin.Copyable
    % Market calculations in Mixed Logit
    %   Used by RCDemand2
    
    properties
        selection
        config
        nind
        iweight
        v
        vx
        x2
        expmu
        edelta % While implementing Hessian at least...
        s
        p
        d
        nonlinprice
        alpha
        sigma
    end
    
    properties
        times_array = @(x, y) bsxfun(@times, x, y);
        plus_array = @(x, y) bsxfun(@plus, x, y);
        transpose_array = @(x)  permute(x, [2,1,3]);
        
        % Transform jacobian to array derivative of vector
        d_array = @(x) reshape(x, [size(x,1), 1, size(x,2)] );
        array2mat = @(x) reshape(x, [size(x,1), size(x,3)] );
    end
    
    methods 
        function [jac, hess] = deltaDerivatives(obj, delta, sigma)
            % Calculates both delta Jacobian and Hessian in BLP
            % Note that matrix inversion is used rather than solving, due to
            % involvement of array mult. As the number of layers is theta, it would be
            % best not to solve.
            [si, dsdd, dsdr] = obj.shareJacobians(delta, sigma);
            [d2sdd2,d2sdr2,d2sddr,d2sdrd] = obj.shareHessians(obj.vx, delta, sigma);
            
            jac = -dsdd\dsdr;
            invdsdd = -inv(dsdd);
            analder = obj.arrayMult(obj.arrayMult(d2sdd2, jac, 2), jac) + ...
                obj.arrayMult(d2sdrd, jac,2) + ...
                obj.arrayMult(d2sddr, jac) + d2sdr2;
            hess = obj.arrayMult(invdsdd, analder);
        end
 
        % Combine vec with inverse share Jacobian.
        % Then multiply with other derivative array to get part of objective
        % Hessian
        function [jac, hess] = objectiveHessian(obj, delta, sigma, vec)
            % Create objective vector Xi
            [~, dsdd, dsdr] = obj.shareJacobians(delta, sigma);
            Xi = dsdd\vec;
            jac = -dsdd\dsdr;
            [d2sdd2,d2sdr2,d2sddr,d2sdrd] = obj.shareHessians(obj.vx, delta, sigma);
            analder = obj.arrayMult(obj.arrayMult(d2sdd2, jac, 2), jac) + ...
                obj.arrayMult(d2sdrd, jac,2) + ...
                obj.arrayMult(d2sddr, jac) + d2sdr2;
            hess = obj.arrayMult(Xi', analder);            
        end
                
        function [si, dsdelta, dssigma] = shareJacobians(obj, delta, sigma)
            % sharecalc can be used instead if nlpart has been invoked
            [~, si] = obj.getShares(delta, sigma);
            wsi = si ./ size(si,2); 
            dssigma = zeros(size(si, 1), length(sigma));
            for k = 1:length(sigma)
                svx = wsi .* obj.vx(:,:,k); % Can be done outside loop
                dssigma(:, k) = sum(svx, 2) - si * sum(svx)';
            end
            dsdelta = diag(sum(wsi, 2)) - wsi*si';
        end

        function [hd, hr, hdr, hrd] = shareHessians(obj, vx, delta, sigma)
            % Calculates the share Hessian arrays as function of delta,
            % arrays of dimensions: MxMxM, MxKxM,MxMxK and MxKxK where
            % M is the number of products in S and K the number of of random coefs.
            % NOTE: The function works, but is slow as it loops over individuals
            % Can thus benefit from C++ calculation. Using symmetry in the outer
            % product also reduces the number of multiplication to a half in 2d and
            % more in 3d.
            [~, si] = obj.getShares(delta, sigma);
            M = size(si, 1);
            V = permute(vx, [1,3,2]);
            
            for i = 1:size(si, 2)
                hess = zeros(M, M, M);
                s = si(:, i);
                ss = s*s';
                own = s - 2 * s .* s; % Alt: own=s with no inner if statement
                for j = 1:M
                    hess(:,:,j) = 2 * ss * s(j);
                    for k = 1:M
                        hess(k, k, j) = hess(k, k, j) - ss(j,k); % deltaJac*s(j);
                        if k ~= j
                            hess(j, k, j) = hess(j, k, j) - ss(j,k);
                            hess(k, j, j) = hess(k, j, j) - ss(j,k);
                        end
                    end
                    hess(j,j,j) = hess(j,j,j) + own(j);
                end
                hessdr = obj.arrayMult(hess, V(:,:,i), 2);
                hessrd = obj.arrayMult(hess, V(:,:,i));
                if i == 1
                    hd = hess; % quad goes here
                    hdr = hessdr; % quad goes here
                    hrd = hessrd; % quad goes here
                    hr = obj.arrayMult(hessdr, V(:,:,i));
                else
                    hd = hd + hess;
                    hdr = hdr + hessdr; % quad goes here
                    hrd = hrd + hessrd; % quad goes here
                    hr = hr + obj.arrayMult(hessdr, V(:,:,i));
                end
            end
            hd = hd/size(si, 2); % Quadrature requires weighting
            hr = hr/size(si, 2);
            hdr = hdr/size(si, 2);
            hrd = hrd/size(si, 2);
        end
    end
    
    methods (Static)        
        function z = arrayMult( x, y, varargin )
            % arrayMult multiply 3d array with matrix/vector
            % Multiplication with vector reduces dimensionality
            % Only two types of post-multiplication have been defined
            % Compare with mmult.m
            isarray = @(a)(length(size(a)) == 3);
            dotmult = nargin == 3 && varargin{1} == 2;
            if isarray(x)
                if dotmult
                    x = permute(x, [1,3,2]);
                end
                sz = size(x);
                sz(2) = size(y, 2);
                z = zeros(sz);
                for i = 1:size(x,3)
                    z(:, :, i) = x(:, :, i) * y;
                end
            elseif isarray(y)
                sz = size(y);
                sz(1) = size(x, 1);
                z = zeros(sz);
                for i = 1:size(y,3)
                    z(:, :, i) = x * y(:, :, i);
                end
            else
                z = x * y;
            end
            if dotmult
                z = permute(z, [1,3,2]);
            end
            if isarray(z)
                if sz(1) == 1
                    z = reshape(z, sz([2,3]) );
                elseif sz(2) == 1
                    z = reshape(z, sz([1,3]) );
                elseif sz(3) == 1
                    z = reshape(z, sz([1,2]) );
                end
            end
        end
        
        % Used in equilibrium finding
        function d = diagArray(x)
            % Create diagonal array from matrix
            % Should perhaps be sparse array
            [m,n] = size(x);
            d = zeros(m,m,n);
            for i = 1:n;
                d(:,:,i) = diag(x(:, i));
            end;
        end
    end
    
    methods
        function obj = RCDemandMarket2(demand)
            % Copy config as a struct to avoid deep copy problems with
            % copy() function. Note that one could keep the common config
            % object, since all instances have same config. 
            obj.config = demand.config.getProperties();
            obj.nind = demand.settings.nind;
            obj.alpha = demand.alpha;
            obj.sigma = demand.sigma;
            obj.config.ces = demand.settings.ces;
            obj.config.compiled = true;
            obj.iweight = demand.iweight;
            obj.nonlinprice = strcmp(demand.var.price, demand.nonlinparams);    
        end
        
        function [s, si] = shares(obj, p)
            if obj.config.ces
                p = log(p);
            end          
            if any(obj.nonlinprice) && any( p ~= obj.p)
                obj.p = p;
                obj.x2(:, obj.nonlinprice) = p;
                % Can update just nonlinprice
                for k = 1:size( obj.x2, 2)
                     obj.vx(:,:,k) = bsxfun(@times, obj.x2(:,k), obj.v(k,:));
                end
            end
            % obj.d set in init() or by user in simulating data:
            edeltapred = exp(obj.d - obj.alpha * p); 
            if obj.config.compiled
                [s, si] = shareCalcCpp(obj.expmu, obj.iweight, edeltapred);
            else
                [s, si] = obj.sharecalc(edeltapred);
            end
        end
       
        function sh = shareJacobian(obj, P)
            if isempty(P) 
                P = obj.p;
            end
            [S, si] = obj.shares(P);
            if any(obj.nonlinprice) 
                theta_p = obj.sigma(obj.nonlinprice);
                % obj.sim.v is the same for all products
                dUi = obj.iweight .* (- obj.alpha + ...
                    theta_p * obj.v(1,:,obj.nonlinprice))';
            else
                dUi = - obj.alpha*obj.iweight;
            end
            sh = zeros(size(si,1));
            for i = 1:size(si,2)
                shi = si(:, i);
                sh = sh + dUi(i) *( diag(shi) - shi*shi' );
            end
            if obj.config.ces
                sh = (sh - diag(S) ) * diag(1./P) ;
            end            
        end
        
        % One can test whether substituting the rhs for all instances of
        % obj.vx(:,:,k) affects performance. 
        function init(obj)
            for k = 1:size(obj.x2, 2)
                obj.vx(:,:,k) = bsxfun(@times, obj.x2(:,k), obj.v(k,:));
            end
            obj.nlpart(obj.sigma);
        end
        
        function der = deltaJacobian(obj, sigma, edelta)
            % Note thate here deltaJacobian takes care of selection
            if obj.config.compiled
                [s, si] = shareCalcCpp(obj.expmu, obj.iweight, edelta(obj.selection));
            else
                [~, si] = obj.sharecalc(edelta(obj.selection));
            end
            wsi =  bsxfun(@times, si , obj.iweight'); % Row means
            dssigma = zeros(size(si, 1), length(sigma));
            for k = 1:length(sigma)
                svx = wsi .* obj.vx(:,:,k); % Can be done outside loop
                dssigma(:, k) = sum(svx, 2) - si * sum(svx)';
            end
            dsdelta = diag(sum(wsi, 2)) - wsi * si';
            % Note that dividing by obj.settings.nind cancels out.
            der = -dsdelta\dssigma;
        end
        
        function [sh, indsh] = getShares(obj, delta, sigma)
            nlpart(obj, sigma); 
            eg = bsxfun(@times, exp(delta), obj.expmu );
            denom = 1 ./ ( 1 + sum(eg));
            indsh = bsxfun(@times, eg , denom);
            sh = indsh * obj.iweight; % Row means
        end
             
        function [sh, indsh] = sharecalc(obj, edelta)
            eg = bsxfun(@times, edelta, obj.expmu );
            denom = 1 ./ ( 1 + sum(eg));
            indsh = bsxfun(@times, eg , denom);
            sh = indsh * obj.iweight; % Row means
        end
        
        function nlpart(obj, sigma)
            [n,m,K] = size(obj.vx);
            mu = zeros(n, m);
            for k = 1:K
                mu = mu + obj.vx(:,:,k) * sigma(k);
            end
            obj.expmu = exp(mu);
        end
        
        function edelta = findDelta(obj, sigma, edelta, tolerance)
            obj.nlpart( sigma );
            if obj.config.compiled
                edelta = findDeltaCpp(obj.s, obj.expmu, obj.iweight, edelta);
            else
                i = 0;
                maxdev = 100;
                while maxdev > tolerance && i < obj.config.fpmaxit
                    eg = bsxfun(@times, edelta, obj.expmu );
                    % Code with intermediate steps
                    % denom = 1 ./ ( 1 + sum(eg));
                    % share = bsxfun(@times, eg, denom) * obj.iweight;
                    % newedel = edelta .* obj.s ./ share;
                    newedel = edelta .* obj.s ./ ( ...
                        bsxfun(@times, eg, 1 ./ ( 1 + sum(eg)) ) ...
                        * obj.iweight);
                    maxdev = max(abs(newedel - edelta));
                    edelta = newedel;
                    i = i + 1;
                end
            end
            obj.edelta = edelta;
        end
          
    end
end

