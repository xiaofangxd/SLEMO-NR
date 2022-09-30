function MOEALPCA(Global)
% <algorithm> <M>
% Multi-objective evolutionary algorithm based on Pareto optimal subspace
% learning

%------------------------------- Reference --------------------------------
% Y. Tian, C. Lu, X. Zhang, K. C. Tan, and Y. Jin, Solving large-scale
% multi-objective optimization problems with sparse optimal solutions via
% unsupervised neural networks, IEEE Transactions on Cybernetics, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Population initialization
    REAL = ~strcmp(Global.encoding,'binary');
    lambdaValue = linspace(0.00001,0.1,Global.N);
    PP = lasso(Global.problem.A,Global.problem.y,'Lambda',lambdaValue,'Alpha',0.5);
    PP = PP';
    PP(PP>=0.5) = 1;
    PP(PP<0.5) = 0;

    m = 10; 
    lamda = 0.1 * norm(Global.problem.A' * Global.problem.y,inf);
%     lamda = 1.0;
    Mask = PP;
    Population = INDIVIDUAL(Mask);
    [Population,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Mask,Global.N);
    
    %% Optimization
    rho = 0.5; 
    symbol = 0;
    while Global.NotTermination(Population)
        Site = rho > rand(1,ceil(Global.N/2));
        if any(Site)
            [lpca,allZero,allOne,symbol] = ModelTraining(Mask(FrontNo==1,:),Mask(FrontNo==2,:),m);
        else
            [lpca,allZero,allOne] = deal([],[],[]);
        end
        MatingPool = TournamentSelection(2,ceil(Global.N/2)*2,FrontNo,-CrowdDis);
        [TotalPopulation,TotalMask] = Operator(Mask(MatingPool,:),lpca,Site,allZero,allOne,symbol,Global.problem.A,Global.problem.y,m,Global.D,Global.N,lamda,Population,Mask);
        
%         Offspring = INDIVIDUAL(OffMask);
        [Population,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(TotalPopulation,TotalMask,Global.N);
%         rho = (rho+sRatio)/2;
        lamda = max(0.98*lamda,1e-2); % (0.85-0.98)
    end
end