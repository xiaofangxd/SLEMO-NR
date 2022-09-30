function [Population,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Mask,N)
% The environmental selection of MOEADSR

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete duplicated solutions
%     success = false(1,length(Population));
    [~,uni] = unique(Population.objs,'rows');
    if length(uni) == 1
        [~,uni] = unique(Population.decs,'rows');
    end
    if length(uni)<N
       add_length = N - length(uni);
       add = randi(length(Population),1,add_length);
       uni = [uni;add'];
    end
    Population = Population(uni);
    Mask       = Mask(uni,:);
    N          = min(N,length(Population));

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;    
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Calculate the ratio of successful offsprings
%     success(uni(Next)) = true;
%     s1     = sum(success(len+1:len+num));
%     s2     = sum(success(len+num+1:end));
%     sRatio = (s1+1)./(s1+s2+1);
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
    Mask       = Mask(Next,:);
end