function [Aa, A] = data_generate_CN(M, W, N, Noise)
Aa = 10.0*rand(M,N);
A = Aa*W;
Aa = Aa + Noise*randn(M,N);
