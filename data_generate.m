function data_generate(M, P_PDG, W, N, r, b, K, Noise, SV)
% play game with other players
s = [0,1;1,0];
ss = zeros(2,2);
% 就这4种情况
% ij策略分别为 00 01 10 11，对应ss(1,1) ss(1,2) ss(2,1) ss(2,2)
for i=1:2
    for j=1:2
        ss(i,j) = s(i,:)*P_PDG*s(:,j);
    end
end


for sv = 1:SV
    sv
    tic
    % Initialize
    for i = 1:N
        node_all(i) = Node_data(N,M,i,sv);
    end
    player = zeros(1,N);
    rand_temp = rand(1,N);
    player(1,rand_temp<=0.5) = 1;
    player(1,rand_temp>0.5) = 0;
    %% caculate payoff for each node
    for t = 1:M
        F = zeros(N,N);
        G = zeros(1,N);
        for i=1:N
            for j=1:N
                F(i,j) = ss(player(i)+1,player(j)+1);
            end
        end
        for i=1:N
            G(1,i) = F(i,:)*W(:,i);
            node_all(i).get_data(G(1,i),sparse(F(i,:)),t);
        end
        
        % update strategies
        temp_player = player;
        for k = 1:N
            [c,~] = find(W(:,k) >= 1);
            if ~isempty(c)
                L = randperm(length(c));
                P = 1.0/(1+exp((G(1,k)-G(1,c(L(1))))/K));
                if rand <= P
                    temp_player(1,k) = player(1,c(L(1)));
                end
            end
        end
        player = temp_player;
    end
    for i=1:N
        save_file = sprintf('.\\data\\EG\\node_%d_sv_%d',i,sv);
        this_node = node_all(i);
        save(save_file,'this_node');
    end
    toc
end