classdef CLPCA < handle
    %Logistic convex PCA
    % X : n*d
    properties(SetAccess = public)
        m = 10;
        k = 20;
        U_save;
    end
    methods
        function obj = CLPCA(m,k)
            obj.m = m;
            obj.k = k;
        end
        
        function train(obj,X)
            [n,d] = size(X);
            u = log((sum(X,1)/n)./(1-sum(X,1)/n));
            Q = 2 * X - 1;
            theta = obj.m * Q;
            L = sum(sum((theta - repmat(u,n,1)).^2));
            [~,~,kV] = svd(theta - repmat(u,n,1));
            U = kV(:,1:obj.k);
            H_1 = U*U';
            H_2 = H_1;
            for t = 1:100
               F = H_1 + (t-2)/ (t+1)*(H_1 - H_2);
               P = 1 ./ (1+exp(-((theta - repmat(u,n,1)) * (U*U') + repmat(u,n,1))));
               temp = 2 * (P - X)'*(theta - repmat(u,n,1));
               del_f = temp + temp' - temp.*eye(size(temp,1));
               tt = F - del_f/L;
               [V,D] = eig(tt);
               lambda = linspace(D(d,d)-obj.k/d,D(1,1),d);
               D(logical(eye(size(D)))) = lambda;
               H_2 = H_1;
               H_1 = V*D*V';
            end
            [V,~] = eig(H_1);
            obj.U_save = V(:,1:obj.k);
        end
        
        function low_dim = reduce(obj,X)
            [n,~] = size(X);
            Q = 2 * X - 1;
            theta = obj.m * Q;
%             p = exp((theta - repmat(obj.u_save,n,1))*obj.U_save)./(1+exp((theta - repmat(obj.u_save,n,1))*obj.U_save));
            p = 1.0 ./ (1+exp(-(theta)*obj.U_save));
            [~,kk] = size(p);
            low_dim = zeros(n,kk);
            rand_temp = rand(n,kk);
            low_dim(rand_temp<=p) = 1;
            low_dim(rand_temp>p) = 0;
        end
        
        function X = recover(obj,low_dim)
            [n,~] = size(low_dim);
            Q = 2 * low_dim - 1;
            theta = obj.m * Q;
%             p = exp(theta * obj.U_save' + repmat(obj.u_save,n,1))./(1+exp(theta * obj.U_save' + repmat(obj.u_save,n,1)));
            p = 1.0 ./ (1+exp(-(theta * obj.U_save')));
            [~,dd] = size(p);
            X = zeros(n,dd);
            rand_temp = rand(n,dd);
            X(rand_temp<=p) = 1;
            X(rand_temp>p) = 0;
        end
    end
end