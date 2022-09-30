% clear all;clc;
function [G] = BA_generate(N, m0, m)
% N=100;%�����ģ
% m0=4;%��ʼ��ȫ��ͨ�����нڵ����
% m=2;%ÿ���һ���ڵ�Ķ�
%%��飬<k>=2*m,һ������m0 >= m

G=zeros(N,N);%�����ɵ�����  
for i=1:m0%���ɳ�ʼȫ��ͨ����
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
    [Ak,index]=sort(Ak);%index�������ԭ����˳����Ϣ, ���൱��ð��         
    Ak_res = zeros(N_now,1);
    for j1=1:N_now
        Ak_res(j1,1)=sum(Ak(1:j1,:));
    end
    for edge=1:m
        pro=rand();
        I=find(Ak_res>=pro);
        posi=index(I(1),1);%�ҵ�Ak_res�е�һ������pro��λ��
        while(G(i+m0,posi) == 1)
            pro=rand();

        I=find(Ak_res>=pro);
        posi=index(I(1),1);
        end
        G(i+m0,posi)=1;
        G(posi,i+m0)=1;            
    end           
end
%G�б�����ǲ���������BA����
%д�ļ�,�����ڽӾ����ʾ
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
    
