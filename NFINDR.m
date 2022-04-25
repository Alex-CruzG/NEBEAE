function [P,indices] = NFINDR(Y,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [P,indices] = NFINDR(Y,N)
%
% N-FINDR endmembers estimation in multi/hyperspectral dataset
%
% Inputs
%   Y --> Multi/hyperspectral dataset as 2D matrix (L x K).
%   N --> Number of endmembers to find.
%
% Outputs
%   P --> Matrix of endmembers (L x N).
%   indices --> Indicies of pure pixels in Y
%
% Bibliographical references:
% [1] Winter, M. E., «N-FINDR: an algorithm for fast autonomous spectral 
%     end-member determination in hyperspectral data», presented at the 
%     Imaging Spectrometry V, Denver, CO, USA, 1999, vol. 3753, págs. 266-275.
%
% DUCD February/2021
% IICO-FC-UASLP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% data size
[L,K] = size(Y);

%% Dimensionality reduction by PCA
U = pca(Y,N);
Yr= U.'*Y;

%% Initialization
Po = zeros(L,N);
indices = zeros(1,K);
IDX = zeros(1,K);
TestMatrix = zeros(N);
TestMatrix(1,:) = 1;
for i = 1:N
    idx = floor(rand*K) + 1;
    TestMatrix(2:N,i) = Yr(1:N-1,idx);
    IDX(i) = idx;
end
actualVolume = abs(det(TestMatrix)); % instead of: volumeactual = abs(det(MatrixTest))/(factorial(p-1));
it = 1;
v1 = -1;
v2 = actualVolume;

%% Algorithm
maxit=3*N;
while it<=maxit && v2>v1
    for k=1:N
        for i=1:K
            actualSample = TestMatrix(2:N,k);
            TestMatrix(2:N,k) = Yr(1:N-1,i);
            volume = abs(det(TestMatrix));  % instead of: volume = abs(det(MatrixTest))/(factorial(p-1));
            if volume > actualVolume
                actualVolume = volume;
                IDX(k) = i;
            else
                TestMatrix(2:N,k) = actualSample;
            end
        end
    end
    it = it+1;
    v1 = v2;
    v2 = actualVolume;
end
for i = 1:N
    Po(:,i) = Y(:,IDX(i));
    indices(IDX(i)) = 1;
end
P=Po./repmat(sum(Po),[L,1]);

function [U] = pca(X, d)
    N = size(X, 2);
    xMean = mean(X, 2);
    XZeroMean = X - repmat(xMean, 1, N);     
    [U,~,~] = svds((XZeroMean*XZeroMean.')/N, d);
return;