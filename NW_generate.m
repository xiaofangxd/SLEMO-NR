%�ȹ���һ����������
%ÿ���ڵ������������ڸ�K/2���ڵ�������KΪż����,���������ӵĸ���Ϊp
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
%�����ڲ����绷����K/2���ڵ��еĵ�һ��
%����֮��ı�
for edge=2:K/2
    for i=1:N-edge
        G(i,i+edge) = 1;%˳ʱ���һ����
        G(i+edge,i) = 1;
        if i > edge%��ʱ���һ����
            G(i-edge,i) = 1;
            G(i,i-edge) = 1;
        else
            G(N+i-edge,i) = 1;
            G(i,N+i-edge) = 1;
        end
    end
end
%�������繹�����
%�������Ը���p��������ѡȡ��NK/2���Խڵ�֮���������ӱ߲����������ǲ������رߺ��Ի�
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
%G�б�����ǲ���������NW����
%д�ļ�,�����ڽӾ����ʾ
% str1 = 'NW_N=';
% str2 = num2str(N);
% %�������������мӱ��ˣ�����˵<k>���ʲô������
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
    
                  
                    

        
    
        