function [P,A,S,Yh]=BEAE12(Yo,n,parameters,Po,oae)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [P,A,S,Yh]=BEAE12(Y,n,parameters,Po,oae)
%
% Estimation by Extended Blind End-member and Abundance Extraction (EBEAE)
% Algorithm and Linear Mixture Model 
%
% D. U. Campos-Delgado et al., "Extended Blind End-Member and Abundance 
%   Extraction for Biomedical Imaging Applications," in IEEE Access, 
%   vol. 7, pp. 178539-178552, 2019, doi: 10.1109/ACCESS.2019.2958985.
%
%
% Input Arguments
%
%   Y = matrix of measurements (MxN)
%   n = order of linear mixture model
%   parameters = 9x1 vector of hyper-parameters in BEAE methodology
%              = [initicond rho lambda epsilon maxiter downsampling  ...
%                      parallel normalization display]
%       initcond = initialization of end-members matrix {1,2,3}
%                                 (1) Maximum cosine difference from mean
%                                      measurement (default)
%                                 (2) Maximum and minimum energy, and
%                                      largest distance from them
%                                 (3) PCA selection + Rectified Linear Unit
%                                 (4) ICA selection (FOBI) + Rectified
%                                 Linear Unit
%                                 (5) N-FINDR endmembers estimation in a 
%                                 multi/hyperspectral dataset (Winter,1999)
%                                 (6) Vertex Component Analysis (VCA)
%                                 (Nascimento and Dias, 2005)
%       rho = regularization weight in end-member estimation 
%             (default rho=0.1);
%       lambda = entropy weight in abundance estimation \in [0,1) 
%                (default lambda=0);
%       epsilon = threshold for convergence in ALS method 
%                 (default epsilon=1e-3); 
%       maxiter = maximum number of iterations in ALS method
%                 (default maxiter=20);
%       downsampling = percentage of reduced random downsampling in end-member 
%                      estimation [0,1) (default downsampling=0.5);
%       parallel = implement parallel computation of abundances (0 -> NO or 1 -> YES)
%                  (default parallel=0);
%       normalization = normalization of estimated end-members (0 -> NO or 1 ->YES)
%                       (default normalization=1);
%       display = show progress of iterative optimization process (0 -> NO or 1 -> YES)
%                 (default display=0);
%   Po = initial end-member matrix (Mxn)
%   oae = only optimal abundance estimation with Po (0 -> NO or 1 -> YES)
%         (default oae = 0)
%
% Output Arguments
%
%   P = matrix of end-members (Mxn)
%   A  = abudances matrix (nxN)
%   S  = scaling vector (Nx1)
%   Yh = estimated matrix of measurements (MxN) --> Yh = P*A*diag(S)
%
% Daniel U. Campos Delgado
% FC-IICO-UASLP
% April/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default hyper-parameters of BEAE algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global NUMERROR

initcond=1;
rho=0.1;
lambda=0;
epsilon=1e-3;
maxiter=20;
downsampling=0.5;
parallel=0;
normalization=1;
display=0;
NUMERROR=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check concistency of input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin~=5
    oae=0;
end
if nargin==0
    disp('The measurement matrix Y has to be used as argument!!');
    return;
elseif nargin==1
    n=2;
end
if nargin==3 || nargin==4 || nargin==5
    if length(parameters)~= 9
        disp('The length of parameters vector is not 9 !!');
        disp('Default values of hyper-parameters are used instead');
    else
        initcond=round(parameters(1));
        rho=parameters(2);
        lambda=parameters(3);
        epsilon=parameters(4);
        maxiter=parameters(5);
        downsampling=parameters(6);
        parallel=parameters(7);
        normalization=parameters(8);
        display=parameters(9);
        if initcond~=1 && initcond~=2 && initcond~=3 && initcond~=4 && initcond~=5 && initcond~=6
            disp('The initialization procedure of end-members matrix is 1,2,3,4,5 or 6!');
            disp('The default value is considered!');
            initcond=1;
        end
        if rho<0
            disp('The regularization weight rho cannot be negative');
            disp('The default value is considered!');
            rho=0.1;
        end
        if lambda<0 || lambda>=1
            disp('The entropy weight lambda is limited to [0,1)');
            disp('The default value is considered!');
            lambda=0;
        end
        if epsilon<0 || epsilon>0.5
            disp('The threshold epsilon cannot be negative or >0.5');
            disp('The default value is considered!');
            epsilon=1e-3;
        end
        if maxiter<0 && maxiter<100
            disp('The upper bound maxiter cannot be negative or >100');
            disp('The default value is considered!');
            maxiter=20;
        end
        if downsampling<0 && downsampling>1
            disp('The downsampling factor cannot be negative or >1');
            disp('The default value is considered!');
            downsampling=0.5;
        end
        if parallel~=0 && parallel~=1
            disp('The parallelization parameter is 0 or 1');
            disp('The default value is considered!');
            parallel=0;
        end
        if normalization~=0 && normalization~=1
            disp('The normalization parameter is 0 or 1');
            disp('The default value is considered!');
            normalization=1;
        end
        if display~=0 && display~=1
            disp('The display parameter is 0 or 1');
            disp('The default value is considered!');
            display=0;
        end
    end
    if n<2
        disp('The order of the linear mixture model has to greater than 2!');
        disp('The default value n=2 is considered!');
        n=2;
    end
end
if nargin==4 || nargin==5
    if ~ismatrix(Po)
        disp('The initial end-members Po must be a matrix !!');
        disp('The initialization is considered by the maximum cosine difference from mean measurement');
        initcond=1;
    else
        if size(Po,1)==size(Yo,1) && size(Po,2)==n
            initcond=0;
        else
            disp('The size of Po must be Mxn!!');
            disp('The initialization is considered based on the input dataset');
            initcond=1;
        end
    end
end
if nargin==5
    if oae~=0 && oae~=1
        disp('The assignment of oae is incorrect!!');
        disp('The initial end-members Po will be improved iteratively from a selected sample');
        oae=0;
    elseif oae==1 && initcond~=0
        disp('The initial end-members Po is not defined properly!');
        disp('Po will be improved iteratively from a selected sample');
        oae=0;
    end
end
if nargin>6
    disp('The number of input arguments is 5 maximum');
    disp('Please check the help documentation');
    return;
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random downsampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismatrix(Yo)
    disp('The measurements matrix Y has to be a matrix');
    return;
end
[M,No]=size(Yo);
if M>No
    disp('The number of spatial measurements has to be larger to the number of time samples!');
    return;
end

I=1:No;
N=round(No*(1-downsampling));
Is=randperm(No,N);
Y=Yo(:,Is);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if normalization==1
    mYm=sum(Y,1);
    mYmo=sum(Yo,1);
else
    mYm=ones(1,N);
    mYmo=ones(1,No);
end
Ym=Y./repmat(mYm,[M 1]);
Ymo=Yo./repmat(mYmo,[M 1]);
NYm=norm(Ym,'fro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection of Initial End-members Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if initcond==1 || initcond==2
    if initcond==1
        Po=zeros(M,n);
        index=1;
        pmax=mean(Yo,2);
        Yt=Yo;
        Po(:,index)=pmax;
    elseif initcond==2
        index=1;
        Y1m=sum(abs(Yo),1);
        [~,Imax]=max(Y1m);
        [~,Imin]=min(Y1m);
        pmax=Yo(:,Imax);
        pmin=Yo(:,Imin);
        K=size(Yo,2);
        II=1:K;
        Yt=Yo(:,setdiff(II,[Imax Imin]));
        Po(:,index)=pmax;
        index=index+1;
        Po(:,index)=pmin;
    end
    while index<n
        ymax=zeros(1,index);
        Imax=zeros(1,index);
        for i=1:index
            e1m=sum(Yt.*repmat(Po(:,i),1,size(Yt,2)),1)./sqrt(sum(Yt.^2,1))./sqrt(sum(Po(:,i).^2,1));
            [ymax(i),Imax(i)]=min(abs(e1m));
        end
        [~,Immax]=min(ymax);
        IImax=Imax(Immax);
        pmax=Yt(:,IImax);
        index=index+1;
        Po(:,index)=pmax;
        II=1:size(Yt,2);
        Yt=Yt(:,setdiff(II,IImax));
    end
elseif initcond==3
    [~,~,VV]=svd(Ym',0);
     W=VV(:,1:n);
     Po=W.*repmat(sign(W'*ones(M,1))',M,1); 
elseif initcond==4
    Yom=mean(Ym,2);
    Yon = Ym - repmat(Yom,1,N);
    [~,S,VV]=svd(Yon',0);
    Yo_w= pinv(sqrtm(S))*VV'*Ym; 
    [V,~,~] = svd((repmat(sum(Yo_w.*Yo_w,1),M,1).*Yo_w)*Yo_w');
    W=VV*sqrtm(S)*V(1:n,:)'; 
    Po=W.*repmat(sign(W'*ones(M,1))',M,1);
elseif initcond==5
    Po=NFINDR(Ym,n);
elseif initcond==6
    Po=VCA(Ym,n);
end      

Po(Po<0)=0;
Po(isnan(Po))=0;
Po(isinf(Po))=0;

if normalization==1
    mPo=sum(Po,1);
    P=Po./repmat(mPo,[M 1]);
else
    P=Po;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternated Least Squares Procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter=1;
J=1e5;
Jp=1e6;
tic;
if display==1
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('EBEAE Linear Unmixing');
        disp(['Model Order =' num2str(n)]);
        if oae==1
            disp('Only the abundances are estimated from Po');
        elseif oae==0 && initcond==0
            disp('The end-members matrix is initialized externally by matrix Po');
        elseif oae==0 && initcond==1
            disp('Po is constructed based on the maximum cosine difference from mean measurement'); 
        elseif oae==0 && initcond==2
            disp('Po is constructed based on the maximum and minimum energy, and largest difference from them');
        elseif oae==0 && initcond==3
            disp('Po is constructed based on the PCA selection + Rectified Linear Unit');
        elseif oae==0 && initcond==4
            disp('Po is constructed based on the ICA selection (FOBI) + Rectified Linear Unit');
        elseif oae==0 && initcond==5
            disp('Po is constructed based on N-FINDR endmembers estimation by Winter (1999)');
         elseif oae==0 && initcond==6
            disp('Po is constructed based on Vertex Component Analysis by Nascimento and Dias (2005)');
        end
end

while (Jp-J)/Jp >= epsilon && iter < maxiter && oae==0 && NUMERROR==0
    
    Am = abundance(Ym,P,lambda,parallel);
    Pp=P;
    if NUMERROR==0
        P = endmember(Ym,Am,rho,normalization); 
    end
    Jp=J;
    J=norm(Ym-P*Am,'fro');
    if J > Jp
        P=Pp; break;
    end
    if display ==1
        disp(['Number of iteration =' num2str(iter)]);
        disp(['Percentage Estimation Error =' num2str(100*J/NYm) '%']);
    end
    iter=iter+1;
    
end
if NUMERROR==0
  
    if oae==1
        A=abundance(Ymo,P,lambda,parallel);        
    else
        Ins=setdiff(I,Is); 
        Ams = abundance(Ymo(:,Ins),P,lambda,parallel);
        A=[Am Ams];
        II=[Is Ins];
        [~,Index]=sort(II);
        A=A(:,Index);
    end
    ElapTime=toc;
    if display ==1
            disp(['Elapsep Time =' num2str(ElapTime)]);
    end
    S=mYmo'; 
    Yh=P*(A.*repmat(mYmo,[n,1]));
else
    disp('Please revise the problem formulation, not reliable results');
    sound(0.1*sin(2*pi*(1:1000)/10))
    P=[];
    S=[];
    A=[];
    Yh=[];
end

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function A = abundance(Y,P,lambda,parallel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A = abundance(Y,P,lambda,parallel)
%
% Estimation of Optimal Abundances in Linear Mixture Model
%
% Input Arguments
% Y = matrix of measurements
% P = matrix of end-members
% lambda =  entropy weight in abundance estimation \in (0,1)
% parallel = implementation in parallel of the estimation
%
% Output Argument
% A = abundances matrix 
%
% Daniel U. Campos-Delgado
% Oct/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global NUMERROR

[M,N]=size(Y);
n=size(P,2);
A=zeros(n,N);

if size(P,1) ~= M
    disp('ERROR: the number of rows in Y and P does not match');
    NUMERROR=1;
    sound(0.1*sin(2*pi*(1:1000)/10))
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute fixed vectors and matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = ones(n,1);
d = 1;
Go=P'*P;
lmin=min(eig(Go));
G=Go-eye(n)*lmin*lambda;
while rcond(G)<1e-6
    lambda=lambda/2;
    G=Go-eye(n)*lmin*lambda;
    if lambda<1e-6
       disp('Unstable numerical results in abundances estimation, update rho!!');
       NUMERROR=1;
       sound(0.1*sin(2*pi*(1:1000)/10))
       return;
    end
end
Gi=eye(n)/G;
T1=Gi*c;
T2=c'*T1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Computation of Abundances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parallel==1
    
    parfor k=1:N

        yk=Y(:,k);
        bk=P'*yk;
        byk=yk'*yk;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute Optimal Unconstrained Solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dk=(bk'*T1-1)/T2;
        ak = Gi*(bk-dk*c);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check for Negative Elements
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if(sum(ak>=0) ~=n)

            Iset = zeros(1,n);

            while(sum(ak<0) ~= 0)    

                Iset(ak<0) = 1;
                L = length(find(Iset));

                Q = n+1+L;
                Gamma = zeros(Q);
                Beta = zeros(Q,1);

                Gamma(1:n,1:n) = G/byk;
                Gamma(1:n,n+1) = c;
                Gamma(n+1,1:n) = c';

                cont = 0;
                for i = 1:n
                    if(Iset(i)~= 0)
                        cont = cont + 1;
                        ind = i; 
                        Gamma(ind,n+1+cont) = 1;
                        Gamma(n+1+cont,ind) = 1;   
                    end
                end

                Beta(1:n) = bk/byk;
                Beta(n+1) = d;
                delta = Gamma\Beta;
                ak = delta(1:n);
                ak(abs(ak)<1e-9) = 0;
            end    
        end
        A(:,k) = single(ak); 
    end
    
else
    
    for k=1:N

        yk=Y(:,k);
        byk=yk'*yk;
        bk=P'*yk;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute Optimal Unconstrained Solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dk=(bk'*T1-1)/T2;
        ak = Gi*(bk-dk*c);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check for Negative Elements
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if(sum(ak>=0) ~=n)

            Iset = zeros(1,n);

            while(sum(ak<0) ~= 0)    

                Iset(ak<0) = 1;
                L = length(find(Iset));

                Q = n+1+L;
                Gamma = zeros(Q);
                Beta = zeros(Q,1);

                Gamma(1:n,1:n) = G/byk;
                Gamma(1:n,n+1) = c;
                Gamma(n+1,1:n) = c';

                cont = 0;
                for i = 1:n
                    if(Iset(i)~= 0)
                        cont = cont + 1;
                        ind = i; 
                        Gamma(ind,n+1+cont) = 1;
                        Gamma(n+1+cont,ind) = 1;   
                    end
                end

                Beta(1:n) = bk/byk;
                Beta(n+1) = d;
                delta = Gamma\Beta;
                ak = delta(1:n);
                ak(abs(ak)<1e-9) = 0;
            end    
        end
        A(:,k) = single(ak); 
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


function P = endmember(Y,A,rho,normalization)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  P = endmember(Y,A,rho,normalization)
%
% Estimation of Optimal End-members in Linear Mixture Model
%
% Input Arguments
% Y = Matrix of measurements
% A =  Matrix of abundances
% rho = Weighting factor of regularization term
% normalization = normalization of estimated profiles (0=NO or 1=YES)
% 
% Output Arguments
% P = Matrix of end-members
%
% Daniel U. Campos-Delgado
% Oct/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global NUMERROR

[n,N]=size(A);
[M,K]=size(Y);
P=zeros(M,n);
R=sum(n - (1:(n-1)));
W=repmat((1./K./sum(Y.^2,1))',1,n);

if size(Y,2) ~= N
    disp('ERROR: the number of columns in Y and A does not match');
    NUMERROR=1;
    sound(0.1*sin(2*pi*(1:1000)/10));
    return;
end

O = single(n*eye(n) - ones(n,n));   
n1 = single(ones(n,1));
m1 = single(ones(M,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Optimal End-members Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T0=A*(W.*A') + rho*O/R;
while rcond(T0)<1e-6
   rho=rho/10;
   T0=A*(W.*A') + rho*O/R;
   if rho<1e-6
       disp('Unstable numerical results in end-members estimation, update rho!!');
       NUMERROR=1;
       return;
   end
end
V = eye(n)/T0;
T2 = Y*(W.*A')*V;
if normalization == 1
    T1 = single(eye(M) - (1/M)*(m1*m1'));
    T3 = (1/M)*m1*n1';
    P_est = T1*T2 + T3;
else
    P_est=T2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate and Project Negative Elements 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_est(P_est<0) = 0;
P_est(isnan(P_est))=0;
P_est(isinf(P_est))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize Optimal Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if normalization==1
    Psum=sum(P_est,1);
    P=P_est./repmat(Psum,M,1);
else 
    P=P_est;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function Po = NFINDR(Y,N)

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
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


function Po = VCA(Y,N)

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
return;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function [U] = pca(X, d)
    N = size(X, 2);
    xMean = mean(X, 2);
    XZeroMean = X - repmat(xMean, 1, N);     
    [U,~,~] = svds((XZeroMean*XZeroMean.')/N, d);
return;
