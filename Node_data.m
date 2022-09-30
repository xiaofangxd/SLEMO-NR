classdef Node_data < handle
    properties(SetAccess = public)
        N = 0;
        M = 0;
%         X = [];
        y = [];
        A = [];
        which_node = 0;
        which_batch = 0;
    end
    
    methods
        function obj = Node_data(N,M,which_node,which_batch)
            obj.N = N;
            obj.M = M;
%             obj.X = X;
            obj.y = zeros(M,1);
            obj.A = sparse(zeros(M,N));
            obj.which_node = which_node;
            obj.which_batch = which_batch;
        end
        
        function get_data(obj,temp_y,temp_A,which_round)
            obj.y(which_round) = temp_y;
            obj.A(which_round,:) = temp_A;
        end
    end

end