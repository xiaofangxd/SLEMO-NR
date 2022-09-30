function [A, y] = data_generate_RN(M, W, N, Noise)
% generate the resistor network
deltaW = 20*rand(1,N);
A = zeros(M,N,N);
for i = 1:N
    for k = 1:M
        for j = 1:N
            A(k,j,i) = sin((1000+deltaW(i))*k)-sin((1000+deltaW(j))*k);
        end
    end

end
for i = 1:N
     y(:,i) = A(:,:,i)*W(:,i);
end
A = A(:,:,:) + Noise*randn(M,N,N);
