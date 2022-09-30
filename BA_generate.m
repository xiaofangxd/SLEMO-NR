% clear all;clc;
function [G] = BA_generate(N, m0, m)
% N=100;%网络规模
% m0=4;%初始的全连通网络中节点个数
% m=2;%每添加一个节点的度
%%这块，<k>=2*m,一般来讲m0 >= m

G=zeros(N,N);%欲生成的网络  
for i=1:m0%生成初始全连通网络
    for j=1:m0
        if i~=j
            G(i,j)=1;
            G(j,i)=1;
        end
    end
end

N_sum = N-m0;
for i=1:N_sum
    N_now=m0+i-1;
    Ak = zeros(N_now,1);
    Ak_tem = zeros(N_now,1);
    degree_now = zeros(N_now,1);
    for j1=1:N_now
         degree_now(j1,1)=sum(G(j1,:));
    end    
    for j=1:N_now
        Ak_tem(j,1) = degree_now(j,1);       
    end
    Ak_sum = sum(Ak_tem,1);
    Ak(:,1)=Ak_tem(:,1)/Ak_sum;
    index=zeros(N_now,1);
    for j=1:N_now
        index(j,1)=j;
    end
    [Ak,index]=sort(Ak);%index保存的是原来的顺序信息, 即相当于冒泡         
    Ak_res = zeros(N_now,1);
    for j1=1:N_now
        Ak_res(j1,1)=sum(Ak(1:j1,:));
    end
    for edge=1:m
        pro=rand();
        I=find(Ak_res>=pro);
        posi=index(I(1),1);%找到Ak_res中第一个大于pro的位置
        while(G(i+m0,posi) == 1)
            pro=rand();

        I=find(Ak_res>=pro);
        posi=index(I(1),1);
        end
        G(i+m0,posi)=1;
        G(posi,i+m0)=1;            
    end           
end
%G中保存的是产生的无向BA网络
%写文件,利用邻接矩阵表示
% str1 = 'BA_N=';
% str2 = num2str(N);
% str3 = '_k=';
% k = m*2;
% str4 = num2str(k);
% str1 = strcat(str1,str2,str3,str4);
% fp = fopen(str1,'w');
% for i=1:N
%     for j=1:N
%         fprintf(fp,'%d',G(i,j));
%     end
%     fprintf(fp,'\n');
% end
% fclose(fp);
% fprintf('Finish of generating BA!\n');
    
