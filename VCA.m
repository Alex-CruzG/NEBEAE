function [P,indices,SNRe] = VCA(Y,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [P,indices,SNRe]=VCA(Y,N)
%
% Vertex Component Analysis algorithm for endmembers estimation in multi/hyperspectral dataset
%  
%
% Inputs
%   Y --> Multi/hyperspectral dataset as 2D matrix (L x K).
%   N --> Number of endmembers to find.
%
% Outputs
%   P --> Matrix of endmembers (L x N).
%   indices --> Indicies of pure pixels in Y
%   SNRe --> SNR estimate of data [dB]
%
% References
%   J. M. P. Nascimento and J. M. B. Dias, ?Vertex component analysis: A 
% fast algorithm to unmix hyperspectral data,? IEEE Transactions on 
% Geoscience and Remote Sensing, vol. 43, no. 4, apr 2005.
%
% DUCD February/2021
% IICO-FC-UASLP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization.
K = size(Y, 2);
L = size(Y, 1);

yMean = mean(Y, 2);
RZeroMean = Y - repmat(yMean, 1, K);
[Ud, ~, ~] = svds(RZeroMean*RZeroMean.'/K, N);
Rd = Ud.'*(RZeroMean);
P_R = sum(Y(:).^2)/K;
P_Rp = sum(Rd(:).^2)/K + yMean.'*yMean;
SNR = abs(10*log10( (P_Rp - (N/L)*P_R) / (P_R - P_Rp) ));
SNRe = SNR;

SNRth = 15 + 10*log(N) + 8;
if (SNR > SNRth) 
    d = N;
    [Ud, ~, ~] = svds((Y*Y.')/K, d);
    Yd = Ud.'*Y;
    u = mean(Yd, 2);
    M =  Yd ./ repmat( sum( Yd .* repmat(u,[1 K]) ) ,[d 1]);
else
    d = N-1;
    r_bar = mean(Y.').';
    Ud = pca(Y, d);
    %Ud = Ud(:, 1:d);
    R_zeroMean = Y - repmat(r_bar, 1, K);
    Yd = Ud.' * R_zeroMean;
     c = zeros(N, 1);
    for j=1:K
        c(j) = norm(Yd(:,j));
    end
    c = repmat(max(c), 1, K);
    M = [Yd; c];
end
e_u = zeros(N, 1);
e_u(N) = 1;
A = zeros(N, N);
% idg - Doesnt match.
A(:, 1) = e_u;
I = eye(N);
k = zeros(K, 1);
for i=1:N
    w = rand(N, 1);
    % idg - Oppurtunity for speed up here.
    tmpNumerator =  (I-A*pinv(A))*w;
    %f = ((I - A*pinv(A))*w) /(norm( tmpNumerator ));
    f = tmpNumerator / norm(tmpNumerator);

    v = f.'*M;
    k = abs(v);
    [~, k] = max(k);
    A(:,i) = M(:,k);
    indices(i) = k;
end
if (SNR > SNRth)
    Po = Ud*Yd(:,indices);
else
    Po = Ud*Yd(:,indices) + repmat(r_bar, 1, N);
end
P=Po./repmat(sum(Po),[L,1]);
return;

function [U] = pca(X, d)
    N = size(X, 2);
    xMean = mean(X, 2);
    XZeroMean = X - repmat(xMean, 1, N);     
    [U,~,~] = svds((XZeroMean*XZeroMean.')/N, d);
return;
