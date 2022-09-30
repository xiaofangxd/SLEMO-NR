function [lpca,allZero,allOne,symbol] = ModelTraining(Mask,add_Mask,m)
% Training RBM and DAE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Determine the size of hidden layers
    if size(Mask,1)==1
        Mask = [Mask;add_Mask];   
    end
    allZero = all(~Mask,1);
    allOne  = all(Mask,1);
    other   = ~allZero & ~allOne;
%     K       = sum(mean(abs(Mask(:,other))>1e-6,1)>rand());
% %     K       = min(max(K,1),size(Mask,1));
%     K       = min(max(K,1));


    %% Train LPCA
%     rbm = RBM(sum(other),K,10,1,0,0.5,0.1);
%     rbm.train(Mask(:,other));
    
    if sum(other)==0
        symbol = 1;% skip reduction
        lpca = LPCA(m);
    else
        symbol = 0;
        lpca = LPCA(m);
        lpca.train(Mask(:,other));
    end
end