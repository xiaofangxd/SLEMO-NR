function [TotalPopulation,TotalMask] = Operator(ParentMask,lpca,Site,allZero,allOne,symbol,A_temp,y_temp,m,D,N,lamda,Population,Mask)
% The operator of MOEA/DSR
%
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Parent1Mask = ParentMask(1:end/2,:);
    Parent2Mask = ParentMask(end/2+1:end,:);
    
    %% Binary variation
    if any(Site)
        if symbol == 1
            OffMask = BinaryVariation(Parent1Mask(Site,:),Parent2Mask(Site,:));
        else
            other   = ~allZero & ~allOne;
            OffTemp = BinaryVariation(lpca.reduce(Parent1Mask(Site,other)),lpca.reduce(Parent2Mask(Site,other)));
            OffTemp = lpca.recover(OffTemp);
            OffMask = false(size(OffTemp,1),size(Parent1Mask,2));
            OffMask(:,other)  = OffTemp;
            OffMask(:,allOne) = true;
        end
    else
        OffMask = [];
    end
    OffMask = [OffMask;BinaryVariation(Parent1Mask(~Site,:),Parent2Mask(~Site,:))];
    
    
    
    %% Local search
    ls = 1;% IST的迭代次数, sometimes 10 is better
    beta = 1;% epsilon, two versions, one is the absolute value, the other is the percentage as in the paper (in line 52)
    T = 5;% neighborhood size
    
    Local_pop = [Population,INDIVIDUAL(OffMask)];   
    [FrontNo,~] = NDSort(Local_pop.objs,Local_pop.cons,length(Local_pop));
    E = Local_pop(FrontNo==1);
    
    for iter = 1:1
    % for iter = 1:0      without local search
    % Generate Archive
    E_obj = cat(1,E.obj);
    f_min = min(E_obj(:,1));
%     f_max = max(temp_obj(:,1));

    % beta = 0.05 * (f_max - f_min);

    E1 = find(E_obj(:,1) < f_min + beta);
%     E2 = find(temp_obj(:,1) >= f_min + beta);
    starting = E1(randi(length(E1)));
    x_k = Local_pop(starting).dec';   
    sparsity = Local_pop(starting).obj(2);
%     if sparsity < 0.2*D
        local_operator = 1;
%     else
%         local_operator = 1;
%     end

    % Select from archive
    
    kl = sparsity - T;
    kr = sparsity + T;
    k__ = sparsity + ceil(T * rand());
    
    % Thresholding search
%     E3 = find(temp_obj(:,2) == k__, 1);
%     if isempty(E3)
%         E3 = find(temp_obj(:,2) >= sparsity & temp_obj(:,2) <= (sparsity+3));
%     end
%     start_next = E3(randi(length(E3)));
    start_next = find(FrontNo(starting)<FrontNo,1);
    if isempty(start_next)
        start_next = randi(length(Local_pop));
    end
    %Local search  
    x_k_1 = Local_pop(start_next).dec';
    
%     if FrontNo(starting) > FrontNo(start_next)
%         [x_k, x_k_1] = deal(x_k_1, x_k);
%     end
    p_k = 1 ./ (1+exp(-(m*(2*x_k-1))));
    p_k_1 = 1 ./ (1+exp(-(m*(2*x_k_1-1))));
    if sum(abs(p_k-p_k_1)) <1e-7
        p_k_1 = p_k_1 + 0.01 * rand(D,1);
    end
    
    for nnn = 1:ls
        if local_operator == 1
            p_k_1(p_k_1==0) = p_k_1(p_k_1==0) + 0.01;
            p_k(p_k==0) = p_k(p_k==0) + 0.01;
            x_k = 0.5 + 0.5/m*log(p_k./(1-p_k));
            x_k_1 = 0.5 + 0.5/m*log(p_k_1./abs(1-p_k_1));
%             alpha_k=((p_k-p_k_1)'*(1*A_temp'*(A_temp*x_k-y_temp) ./ (p_k - p_k.^2)-1*A_temp'*(A_temp*x_k_1-y_temp) ./ (p_k_1 - p_k_1.^2)))/((p_k-p_k_1)'*(p_k-p_k_1));           
%             alpha_k = 2;
            s_t = 0.5/m*(log((p_k .* (1-p_k_1) ./ (1-p_k) ./ (p_k_1))));
            r_t = (1*A_temp'*(A_temp*x_k-y_temp) ./ (p_k - p_k.^2)) - (1*A_temp'*(A_temp*x_k_1-y_temp) ./ (p_k_1 - p_k_1.^2));
            alpha_k = (s_t' * r_t) / (s_t' * s_t);
            
            if alpha_k < 0
                ccc = 1;
            end
            
            if alpha_k == 0
                alpha_k = alpha_k + 1e-4;
            end
            
%             lamda = 1.0*rand();
% lamda = 0;
            
            u_k=p_k-(1/alpha_k)*(1*A_temp'*(A_temp*x_k-y_temp) ./ (p_k - p_k.^2))/2/m;
        
            findbigger=(abs(u_k)-lamda/alpha_k)>0;
            new_p_k = zeros(1,size(x_k,1));
            new_p_k(findbigger)=sign(u_k(findbigger)).*(abs(u_k(findbigger))-lamda/alpha_k);
%             findsmaller=(abs(u_k)-lamda/alpha_k)<=0;
%             new_p_k(findsmaller)=0;
            
            p_k_1 = p_k;
            p_k = new_p_k';
            if sum(abs(p_k-p_k_1)) <1e-4
                p_k_1 = p_k_1 + 0.01 * rand(D,1);
            end
        elseif local_operator == 2
            % elastic net
            alpha = 0.001;
            rho = 1/3;
            a = sum(A_temp.^2,1)/length(y_temp) + alpha*(1-rho);
            b = sum((A_temp.^2).*p_k',1)/length(y_temp);
            
            u_k = (b./a);
            alpha_k = (alpha*rho./a);
            findbigger=(abs(b)-alpha*rho)>0;
            new_p_k = zeros(1,D);
            new_p_k(findbigger)=sign(u_k(findbigger)).*(abs(u_k(findbigger))-alpha_k(findbigger));
            findsmaller=(abs(b)-alpha*rho)<=0;
            new_p_k(findsmaller)=0;
            
            p_k = new_p_k';
        elseif local_operator == 3
            % ITH/L0.5
            x_k = 0.5 + 0.5/m*log(p_k./(1-p_k));
            
            p_k = p_k - 1*A_temp'*(A_temp*x_k-y_temp) ./ (p_k-p_k.^2);
            [sorted_p_k,~] = sort(abs(p_k),'descend');
            lambda = sorted_p_k(k__+1)^1.5 * 96^0.5/9;
            findbigger = abs(p_k) > 54^1.5/4*lambda^(2/3);
            new_p_k = zeros(1,D);
            phi = 2/3*p_k.*(1+cos((2*3.14/3)-(2/3*acos(lambda/8*(abs(p_k)/3).^(-1.5)))));
            new_p_k(findbigger) = phi(findbigger);
            
            p_k = new_p_k';
        else
            new_p_k = lasso(A_temp,y_temp,'Lambda',0.02,'Alpha',0.5);
        end
    end
    % Truncation by multilevel
    E_temp = [];
    z = new_p_k;
    E_temp = [E_temp;z];
%     c = 1;
%     k__ = sum(z~=0);
    if sum(z~=0) > kr
       [~,I] = sort(z);
       z(I(1:sum(z~=0)-kr)) = 0;
    end
    while sum(z~=0) >= kl
        if sum(z~=0) == 0
            break;
        end
        z(z == min(z(z~=0))) = 0;%每次把非0值当中最小的置0
        E_temp = [E_temp;z];
%         c = c + 1;
    end
    New_P = E_temp;
    % 归一化
    max_P = max(New_P,[],2);
    min_P = min(New_P,[],2);
    New_P = (New_P - repmat(min_P,1,D)) ./ (max_P - min_P + 1e-4);
    % Resample
    local_pop_size = size(New_P,1);
    rnd_tmp = rand(local_pop_size,D);
    LocalMask = zeros(local_pop_size,D);
    LocalMask(rnd_tmp<New_P) = 1;
    
    % Update Archive  
    New_E = INDIVIDUAL(LocalMask);
    New_E_obj = cat(1,New_E.obj);
    
    % step 1: update the external archive    
    E_bigger = New_E(New_E_obj(:,1)>=f_min+beta);
    E_smaller = New_E(New_E_obj(:,1)<f_min+beta);
% 1.1 for the smaller
    if ~isempty(E_smaller)
        tmp1 = cat(1,E_smaller.obj);
        [~,ia,ib] = intersect(E_obj(:,2),tmp1(:,2));
        if ~isempty(ia)
            % for the overlapped, replace them if with smaller f1
            smaller = E_obj(ia,1)>tmp1(ib,1);
            E(ia(smaller)) = E_smaller(ib(smaller));
            % for the not overlapped, add to E directly
            E_smaller(ib) = [];
            E = [E,E_smaller];
        else
            E = [E,E_smaller];
        end
    end
%1.2 for the bigger
    if ~isempty(E_bigger)
        % add to E and remove dominated
        E = [E,E_bigger];
        [FrontNoE,~] = NDSort(E.objs,E.cons,length(E));
        E = E(FrontNoE == 1);
    end

    % step 2: reduce the size of the external archive
    E_obj = cat(1,E.obj);%此时的E已经变了
    f_min = min(E_obj(:,1));
%     f_max = max(temp_obj(:,1));

    E1 = E(E_obj(:,1) < f_min + beta);
    E2 = E(E_obj(:,1) >= f_min + beta);
    if length(E1) > T
       E1_obj = cat(1,E1.obj);
       E1_obj = E1_obj(:,2);
       [~,I] = sort(E1_obj);
       E1(I(T:end)) = [];
       beta = max(0.8*beta,1e-6);
    end
    if length(E2) > T
        E2_obj = cat(1,E2.obj);
        E2_obj = E2_obj(:,2);
        [~,I] = sort(E2_obj);
        E2(I(1:length(E2)-T)) = [];
    end
    E = [E1,E2];
    end
   
   TotalPopulation = [Local_pop,E];
%     TotalPopulation = [Local_pop,New_E];
    TotalMask = cat(1,TotalPopulation.dec);

    
end

function Offspring = BinaryVariation(Parent1,Parent2)
% One point crossover and bitwise mutation

    [proC,proM] = deal(1,1);
    [N,D] = size(Parent1);
    k = repmat(1:D,N,1) > repmat(randi(D,N,1),1,D);
    k(repmat(rand(N,1)>proC,1,D)) = false;
    Offspring1    = Parent1;
    Offspring2    = Parent2;
    Offspring1(k) = Parent2(k);
    Offspring2(k) = Parent1(k);
    
    Site1 = rand(N,D) < proM/D;
    Offspring1(Site1) = ~Offspring1(Site1);
    Site2 = rand(N,D) < proM/D;
    Offspring2(Site2) = ~Offspring2(Site2);
    
    % variable swap
%     swap_indicator = (rand(N,D) >= 0.5);
%     temp = Offspring2(swap_indicator);
%     Offspring2(swap_indicator) = Offspring1(swap_indicator);
%     Offspring1(swap_indicator) = temp;
    
    Offspring     = [Offspring1;Offspring2];
    
end

