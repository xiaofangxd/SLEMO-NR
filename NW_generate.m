%先构造一个规则网络
%每个节点与它左右相邻各K/2个节点相连（K为偶数）,边重新连接的概率为p
function [G] = NW_generate(N, K, p)
% K = 4;
% N = 100;
% p = 0.1;
G = zeros(N,N);
for i=1:N-1
    G(i,i+1) = 1;
    G(i+1,i) = 1;
end
G(1,N) = 1;G(N,1) = 1;
%构造内层网络环，即K/2个节点中的第一个
%构造之后的边
for edge=2:K/2
    for i=1:N-edge
        G(i,i+edge) = 1;%顺时针的一条边
        G(i+edge,i) = 1;
        if i > edge%逆时针的一条边
            G(i-edge,i) = 1;
            G(i,i-edge) = 1;
        else
            G(N+i-edge,i) = 1;
            G(i,N+i-edge) = 1;
        end
    end
end
%规则网络构造完成
%接下来以概率p在网络中选取（NK/2）对节点之间进行随机加边操作，条件是不产生重边和自环
for i=1:N*K/2
    node1 = round(rand()*N);
    while(~node1)
        node1 = round(rand()*N);
    end
    node2 = round(rand()*N);
    while(~node2 || node1 == node2 || G(node1,node2))
        node2 = round(rand()*N);
    end
    if(rand() < p)
        G(node1,node2) = 1;G(node2,node1) = 1;
    end
end    
%G中保存的是产生的无向NW网络
%写文件,利用邻接矩阵表示
% str1 = 'NW_N=';
% str2 = num2str(N);
% %这里面向网络中加边了，不好说<k>变成什么样子了
% % str3 = '_k=';
% % str4 = num2str(K);
% str1 = strcat(str1,str2);
% fp = fopen(str1,'w');
% for i=1:N
%     for j=1:N
%         fprintf(fp,'%d',G(i,j));
%     end
%     fprintf(fp,'\n');
% end
% fclose(fp);
% fprintf('Finish of generating NW!\n');
    
                  
                    

        
    
        