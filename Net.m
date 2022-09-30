cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));
% Generate time series for complex networks
clc,clear;
% rand('twister',sum(100*clock));
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));

SV = 4;
M = 10;
batch_size = 24; % ÿ�ζ�ȡbatch_size���ڵ������
% Node = 5000;
generate_new_data = 0;
which_algorithm = 1;% 1 for mine and 2 for others
sizes = 0;
% load('W.mat');
% %% Load mtx

% data = importdata('USAir97.mtx'); % 332   
% data = importdata('Harvard500.mtx'); % 500
% data = importdata('celegansneural.mtx');
% data = importdata('130bit.mtx'); % 584
% data = importdata('bibd_12_5.mtx'); % 792
% data = importdata('lp_beaconfd.mtx'); % 295
 data = importdata('polbooks.mtx'); % 105
% data = importdata('TSOPF.mtx'); % 5374


Node = data.data(1,1);
data.data(1,:) = [];
W = zeros(Node,Node);
for i = 1:size(data.data,1)
    W(data.data(i,1),data.data(i,2)) = 1;
%     W(data.data(i,2),data.data(i,1)) = 1;
end
for i=1:Node
   W(i,i) = 0; 
end
clear data;

% Load edges
% data = importdata('infect-dublin.edges');%node 410
% data = importdata('email-Eu-core.txt');%node 1005 edges 25571
% Node = 1005;
% W = zeros(Node,Node);
% for i = 1:size(data,1)
%     W(data(i,1)+1,data(i,2)+1) = 1;
% %     W(data(i,2),data(i,1)) = 1;
% end
% for i = 1:Node
%     W(i,i)=0;
% end
% clear data;


%% parameter setting

r = 0.7;b = 1.2;
P_PDG = [1 0;b 0];
P_SG = [1 1-r;1+r 0];
K = 0.1;

Noise = 0.00;
pc = 0.06;
pop_size = 100;
re_W = [];
pop_all = [];
fit_all = [];

% generate complex networks (BA, ER, WS)
if sizes == 1
    W = ER_generate(Node, 0.1);
elseif sizes == 2
     W = BA_generate(Node, 25, 25);
elseif sizes == 3
     W = NW_generate(Node, 35, 0.5);
elseif sizes == 4
    W = WS_generate(Node, 50, 1);
end
% W = W.*(1.0+(6.0-1.0).*rand(N,N));
% initialize the value of Nodes
% generate training data
y = [];
AA = [];
model_type = 0;
%% Generate data and save
if generate_new_data == 1

    data_generate(M, P_PDG, W, Node, r, b, K, Noise, SV);
%     [A, Aa] = data_generate_RN(M, W, Node, Noise);
%     [A, Aa] = data_generate_CN(M, W, Node, Noise);model_type = 3;
%     y = [y;Aa];
%     AA = [AA;A];

% % �������������
% for i = 1:Node
%     if mod(i,100)==1
%         batch = ceil(i/100);
%         savefile=sprintf('N_%d_batch_%d_saveddata',Node,batch);
%         if batch*100 <= Node
%             temp_A = AA(:,:,(1+100*(batch-1)):100*batch);
%             temp_y = y(:,(1+100*(batch-1)):100*batch);
%             save(savefile,'temp_A','temp_y');
%         else
%             temp_A = AA(:,:,(1+100*(batch-1)):end);
%             temp_y = y(:,(1+100*(batch-1)):end);
%             save(savefile,'temp_A','temp_y');
%         end
%     end
% end
% clear AA y A Aa;
end


% Continue
% Should also change the following bt
% bt = 111:ceil(Node/batch_size);
% contifile = sprintf('.\\save_W\\temp\\re_W_83.mat');
% temp_W = load(contifile);
% re_W = temp_W.re_W;


tic
%% Algorithm  starts here
tic
for bt = 1:ceil(Node/batch_size)
    bat_start = (bt-1)*batch_size+1
    bat_end = bt*batch_size;
    if bat_end > Node
        bat_end = Node;
    end
    this_batch_size = bat_end-bat_start+1;
    y = zeros(SV*M,this_batch_size);
    A = zeros(SV*M,Node,this_batch_size);
    for i = bat_start:bat_end
        for j = 1:SV
            loadfile=sprintf('.\\data\\EG\\node_%d_sv_%d',i,j);
            load(loadfile);
            y(((j-1)*M+1):(j*M),mod(i-1,batch_size)+1) = this_node.y;
            A(((j-1)*M+1):(j*M),:,mod(i-1,batch_size)+1) = full(this_node.A);
        end
    end
    
    parfor nn = 1:this_batch_size
            Global = GLOBAL('-algorithm',@MOEALPCA,'-problem',{@NetworkReconstruction,A(:,:,nn),y(:,nn)},'-N',pop_size,'-evaluation',2e4,'-D',Node);
            Global.Start();
            result = Global.result(end,end);
            result = result{1};
            re_W = [re_W;result(1,1).dec];
%             pop_all = [pop_all;result];
            
        end 
    end
    savefile = sprintf('.\\save_W\\re_W_%d',bt);
    save(savefile,'re_W');
end
toc
re_W = re_W';
for i=1:Node
   re_W(i,i) = 0; 
end
% fitness = zeros(Node, pop_size, 2);
% for i=1:Node
%    fitness(i,:,:) = fit_all(:, 2*i-1:2*i);
% end
toc

save('W.mat','W');
save('re_W.mat','re_W');
% save('non_dom','non_dom');
% save('fitness', 'fitness');
% save('fit_all', 'fit_all');
% if which_algorithm == 1
%     save('pop_all_7.mat','pop_all');
% elseif which_algorithm == 2
%     save('fit_all', 'fit_all');
end