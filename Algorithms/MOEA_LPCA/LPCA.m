classdef LPCA < handle
    %Logistic PCA
    % X : n*d
    properties(SetAccess = public)
        m = 10;
        k = 20;
        u_save;
        U_save;
    end
    methods
        function obj = LPCA(m)
            obj.m = m;
        end
        
        function train(obj,X)
            [n,d] = size(X);
            Q = 2 * X - 1;
            %             u = log((sum(X,1)/n)./(1-sum(X,1)/n));
            u = sum(obj.m*Q,1);
            
            %             [~,~,kV] = svd(Q);
            [~,~,kV] = svd(obj.m * Q - repmat(u,n,1));
            %             obj.k = min(n,d);
            if n < d
                obj.k = n;
                U = kV(:,1:obj.k);
                for iter = 1:10
                    theta = (obj.m * Q - repmat(u,n,1)) * (U * U') + repmat(u,n,1);
                    Z = theta + 4 * (X - 1./(1+exp(-theta)));
                    u = sum(Z - obj.m * Q * (U * U'),1)/n;
                    %                 tt = (obj.m * Q - repmat(u,n,1))' * (Z - repmat(u,n,1)) + (Z - repmat(u,n,1))' * (obj.m * Q - repmat(u,n,1)) - (obj.m * Q - repmat(u,n,1))'*(obj.m * Q - repmat(u,n,1));
                    [E,~] = eig((obj.m * Q - repmat(u,n,1))' * (Z - repmat(u,n,1)) + (Z - repmat(u,n,1))' * (obj.m * Q - repmat(u,n,1)) - (obj.m * Q - repmat(u,n,1))'*(obj.m * Q - repmat(u,n,1)));
                    U = E(:,1:obj.k);
                end
                obj.u_save = u;
                obj.U_save = U;
            else                
                count = 0;
                for temp_k = d:-1:1
                    obj.k = temp_k;
                    U = kV(:,1:obj.k);
                    
                    for iter = 1:10
                        theta = (obj.m * Q - repmat(u,n,1)) * (U * U') + repmat(u,n,1);
                        Z = theta + 4 * (X - 1./(1+exp(-theta)));
                        u = sum(Z - obj.m * Q * (U * U'),1)/n;
                        %                 tt = (obj.m * Q - repmat(u,n,1))' * (Z - repmat(u,n,1)) + (Z - repmat(u,n,1))' * (obj.m * Q - repmat(u,n,1)) - (obj.m * Q - repmat(u,n,1))'*(obj.m * Q - repmat(u,n,1));
                        [E,~] = eig((obj.m * Q - repmat(u,n,1))' * (Z - repmat(u,n,1)) + (Z - repmat(u,n,1))' * (obj.m * Q - repmat(u,n,1)) - (obj.m * Q - repmat(u,n,1))'*(obj.m * Q - repmat(u,n,1)));
                        U = E(:,1:obj.k);
                    end
                    
                    theta = (obj.m * Q - repmat(u,n,1)) * (U * U') + repmat(u,n,1);
                    D1 = sum(sum(-2 * X.* theta + 2 * log(1 + exp(theta))));
                    D2 = sum(sum(-2 * X.* repmat(u,n,1) + 2 * log(1 + exp(repmat(u,n,1)))));
                    if (1-D1/D2) <= 0.95 || count == 10% 如果不满足，继承上一代的结果并跳出
                        if temp_k == d % 不降维
                            break;
                        end
                        u = temp_u;
                        U = temp_U;
                        break;
                    else %如果满足，保存这一代的结果
                        count = count + 1;
                        temp_k = ceil(temp_k/2);
                        temp_u = u;
                        temp_U = U;
                    end
                    
                end
                obj.u_save = u;
                obj.U_save = U;
            end
        end
        function low_dim = reduce(obj,X)
            [n,~] = size(X);
            Q = 2 * X - 1;
            theta = obj.m * Q;
%             p = exp((theta - repmat(obj.u_save,n,1))*obj.U_save)./(1+exp((theta - repmat(obj.u_save,n,1))*obj.U_save));
            p = 1.0 ./ (1+exp(-(theta - repmat(obj.u_save,n,1))*obj.U_save));
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
            p = 1.0 ./ (1+exp(-(theta * obj.U_save' + repmat(obj.u_save,n,1))));
            [~,dd] = size(p);
            X = zeros(n,dd);
            rand_temp = rand(n,dd);
            X(rand_temp<=p) = 1;
            X(rand_temp>p) = 0;
        end
    end
end