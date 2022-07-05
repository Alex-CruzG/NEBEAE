function [P,A,Ds,S,Yh]=NEBEAE(Yo,N,parameters,Po,oae)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [P,A,Ds,S,Yh]=NEBEAE(Y,N,parameters,Po,oae)
%
% Estimation of Optimal End-members and Abundances in Nonlinear Mixture
%  Model (Normalized Model)
%
%
% Input Arguments
%
%   Y = matrix of measurements (LxK)
%   n = order of linear mixture model
%   parameters = 8x1 vector of hyper-parameters in BEAE methodology
%              = [initicond rho lambda epsilon maxiter downsampling  ...
%                      parallel display]
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
%       downsampling = percentage of random downsampling in end-member 
%                      estimation [0,1) (default downsampling=0.5);
%       parallel = implement parallel computation of abundances (0 -> NO or 1 -> YES)
%                  (default parallel=0);
%       display = show progress of iterative optimization process (0 -> NO or 1 -> YES)
%                 (default display=0);
%   Po = initial end-member matrix (LxN)
%   oae = only optimal abundance estimation with Po (0 -> NO or 1 -> YES)
%         (default oae = 0)
%
% Output Arguments
%
%   P = matrix of end-members (LxN)
%   A  = abudances matrix (NxK)
%   Ds = vector of nonlinear interaction levels (Kx1)
%   S  = scaling vector (Kx1)
%   Yh = estimated matrix of measurements (LxK) --> Yh = P*A*diag(S)
%
% Daniel U. Campos Delgado
% August/2021
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
    N=2;
end
if nargin==3 || nargin==4 || nargin==5
    if length(parameters)~= 8
        disp('The length of parameters vector is not 8 !!');
        disp('Default values of hyper-parameters are used instead');
    else
        initcond=round(parameters(1));
        rho=parameters(2);
        lambda=parameters(3);
        epsilon=parameters(4);
        maxiter=parameters(5);
        downsampling=parameters(6);
        parallel=parameters(7);
        display=parameters(8);
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
        if display~=0 && display~=1
            disp('The display parameter is 0 or 1');
            disp('The default value is considered!');
            display=0;
        end
    end
    if N<2
        disp('The order of the linear mixture model has to greater than 2!');
        disp('The default value n=2 is considered!');
        N=2;
    end
end
if nargin==4 || nargin==5
    if ~ismatrix(Po)
        disp('The initial end-members Po must be a matrix !!');
        disp('The initialization is considered by the maximum cosine difference from mean measurement');
        initcond=1;
    else
        if size(Po,1)==size(Yo,1) && size(Po,2)==N
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
[L,Ko]=size(Yo);
if L>Ko
    disp('The number of spatial measurements has to be larger to the number of time samples!');
    return;
end

I=1:Ko;
K=round(Ko*(1-downsampling));
Is=randperm(Ko,K);
Y=Yo(:,Is);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mYm=sum(Y,1);
mYmo=sum(Yo,1);
Ym=Y./repmat(mYm,[L 1]);
Ymo=Yo./repmat(mYmo,[L 1]);
NYm=norm(Ym,'fro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection of Initial End-members Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if initcond==1 || initcond==2
    if initcond==1
        Po=zeros(L,N);
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
    while index<N
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
     W=VV(:,1:N);
     Po=W.*repmat(sign(W'*ones(L,1))',L,1); 
elseif initcond==4
    Yom=mean(Ym,2);
    Yon = Ym - repmat(Yom,1,N);
    [~,S,VV]=svd(Yon',0);
    Yo_w= pinv(sqrtm(S))*VV'*Ym; 
    [V,~,~] = svd((repmat(sum(Yo_w.*Yo_w,1),M,1).*Yo_w)*Yo_w');
    W=VV*sqrtm(S)*V(1:N,:)'; 
    Po=W.*repmat(sign(W'*ones(L,1))',L,1);
elseif initcond==5
    Po=NFINDR(Ym,N);
elseif initcond==6
    Po=VCA(Ym,N);
end      


Po(Po<0)=0;
Po(isnan(Po))=0;
Po(isinf(Po))=0;
mPo=sum(Po,1);
P=Po./repmat(mPo,[L 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternated Least Squares Procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter=1;
J=1e5;
Jp=1e6;

Dm=0*ones(K,1);

if display==1
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('EBEAE Linear Unmixing');
        disp(['Model Order =' num2str(N)]);
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
    
    Am = abundance(Y,Ym,P,Dm,lambda,parallel);
    Dm = probanonlinear(Y,Ym,P,Am,parallel);
    Pp = P;
    if NUMERROR==0
        P = endmember(Y,Ym,Pp,Am,Dm,rho); 
    end

    Jp=J;
    J=norm(Ym - repmat((1-Dm)',[L,1]).*(P*Am) - repmat(Dm',[L,1]).*((P*Am).*Y),'fro');
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
        J=1e5;
        Jp=1e6;
        D=zeros(K,1);
        iter=1;
        
        while (Jp-J)/Jp >= epsilon && iter < maxiter
            
            A=abundance(Yo,Ymo,P,D,lambda,parallel);
            D=probanonlinear(Yo,Ymo,P,A,parallel);
            
            Jp=J;
            J=norm(Ymo - repmat((1-D)',[L,1]).*(P*A) - repmat(D',[L,1]).*((P*A).*Yo),'fro');
            iter=iter+1;
            
        end
        disp(['Percentage Estimation Error =' num2str(100*J/NYm) '%']);
    else
        Ins=setdiff(I,Is); 
        J=1e5;
        Jp=1e6;
        Dms=zeros(length(Ins),1);
        iter=1;
        Ymos=Ymo(:,Ins);
        Yos=Yo(:,Ins);
       
        while (Jp-J)/Jp >= epsilon && iter < maxiter
            
            Ams=abundance(Yos,Ymos,P,Dms,lambda,parallel);
            Dms=probanonlinear(Yos,Ymos,P,Ams,parallel);
            Jp=J;
            J=norm(Ymos - repmat((1-Dms)',[L,1]).*(P*Ams) - repmat(Dms',[L,1]).*((P*Ams).*Yos),'fro');
            iter=iter+1;
            
        end        
        
        A=[Am Ams];
        D=[Dm; Dms];
        II=[Is Ins];
        [~,Index]=sort(II);
        A=A(:,Index);
        D=D(Index);
    end
    ElapTime=toc;
    if display ==1
            disp(['Elapsep Time =' num2str(ElapTime)]);
    end
    S=mYmo';
    AA=(A.*repmat(mYmo,[N,1]));
    Yh=repmat((1-D)',[L,1]).*(P*AA) + repmat(D',[L,1]).*((P*AA).*Yo);
    Ds=probanonlinear(Yo,Yo,P,A,parallel);
else
    disp('Please revise the problem formulation, not reliable results');
    sound(0.1*sin(2*pi*(1:1000)/10))
    P=[];
    Ds=[];
    S=[];
    A=[];
    Yh=[];
end

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function A = abundance(Z,Y,P,D,lambda,parallel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A = abundance(Z,Y,P,D,lambda,parallel)
%
% Estimation of Optimal Abundances in Nonlinear Mixture Model
%
% Input Arguments
% Z --> matrix of measurements
% Y --> matrix of normalized measurements
% P --> matrix of end-members
% D --> vector of probabilities of nonlinear mixing
% lambda -->  entropy weight in abundance estimation \in (0,1)
% parallel --> implementation in parallel of the estimation
%
% Output Argument
% A = abundances matrix 
%
% Daniel U. Campos-Delgado
% February/2021
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Computation of Abundances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parallel==1
    
    parfor k=1:N

        yk=Y(:,k);
        zk=Z(:,k);
        byk=yk'*yk;
        dk=D(k);
        deltakn=(1-dk)*ones(n,1)+dk*P'*zk;
        Pk=P.*((1-dk)*ones(M,n)+dk*zk*ones(1,n));
        bk=Pk'*yk;
        Go=Pk'*Pk;
        eGo=eig(Go);
        eGo(isnan(eGo))=1e6;
        eGo(isinf(eGo))=1e6;
        lmin=min(eGo);
        G=Go-eye(n)*lmin*lambda;
        Gi=eye(n)/G;
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute Optimal Unconstrained Solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sigma=(deltakn'*Gi*bk-1)/(deltakn'*Gi*deltakn);
        ak = Gi*(bk-deltakn*sigma);

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

                Gamma(1:n,1:n) = G;
                Gamma(1:n,n+1) = deltakn*byk;
                Gamma(n+1,1:n) = deltakn';

                cont = 0;
                for i = 1:n
                    if(Iset(i)~= 0)
                        cont = cont + 1;
                        ind = i; 
                        Gamma(ind,n+1+cont) = 1;
                        Gamma(n+1+cont,ind) = 1;   
                    end
                end

                Beta(1:n) = bk;
                Beta(n+1) = 1;
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
        zk=Z(:,k);
        byk=yk'*yk;
        dk=D(k);
        deltakn=(1-dk)*ones(n,1)+dk*P'*zk;
        Pk=P.*((1-dk)*ones(M,n)+dk*zk*ones(1,n));
        bk=Pk'*yk;
        Go=Pk'*Pk;
        eGo=eig(Go);
        eGo(isnan(eGo))=1e6;
        eGo(isinf(eGo))=1e6;
        lmin=min(eGo);
        G=Go-eye(n)*lmin*lambda;
        Gi=eye(n)/G;
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute Optimal Unconstrained Solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sigma=(deltakn'*Gi*bk-1)/(deltakn'*Gi*deltakn);
        ak = Gi*(bk-deltakn*sigma);

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

                Gamma(1:n,1:n) = G;
                Gamma(1:n,n+1) = deltakn*byk;
                Gamma(n+1,1:n) = deltakn'';

                cont = 0;
                for i = 1:n
                    if(Iset(i)~= 0)
                        cont = cont + 1;
                        ind = i; 
                        Gamma(ind,n+1+cont) = 1;
                        Gamma(n+1+cont,ind) = 1;   
                    end
                end

                Beta(1:n) = bk;
                Beta(n+1) = 1;
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

function D=probanonlinear(Z,Y,P,A,parallel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  P = probanonlinear(Z,Y,P,A,parallel)
%
% Estimation of Probability of Nonlinear Mixtures 
%
% Input Arguments
% Z --> matrix of measurements
% Y -> matrix of normalized measurements
% P --> matrix of end-members
% A -->  matrix of abundances
% parallel = implementation in parallel of the estimation
% 
% Output Arguments
% D = Vector of probabilities of Nonlinear Mixtures
%
% Daniel U. Campos-Delgado
% February/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=size(Y,2);
D=zeros(K,1);
G=sum(Z);

if parallel==1
    
    parfor k=1:K

        yk=Y(:,k);
        zk=Z(:,k);
        ak=A(:,k);
        ek=P*ak;        
        T1=ek-yk;
        T2=ek-ek.*zk;
        dk=min([1 T1'*T2/(T2'*T2)]);
        D(k)=dk;
        
    end
    
else
    
    for k=1:K

        yk=Y(:,k);
        zk=Z(:,k);
        ak=A(:,k);
        ek=P*ak;   
        T1=ek-yk;
        T2=ek-ek.*zk;
        dk=min([1 T1'*T2/(T2'*T2)]);
        D(k)=dk;
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


function P = endmember(Z,Y,Po,A,D,rho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  P = endmember(Z,Y,P,A,rho,normalization)
%
% Estimation of Optimal End-members in Linear Mixture Model
%
% Input Arguments
% Z --> matrix of measurements
% Y -> matrix of normalized measurements
% P --> matrix of end-members
% A -->  matrix of abundances
% D --> vector of probabilities of nonlinear mixture
% rho = Weighting factor of regularization term
% normalization = normalization of estimated profiles (0=NO or 1=YES)
% 
% Output Arguments
% P --> matrix of end-members
%
% Daniel U. Campos-Delgado
% Feb/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Gradient of Cost Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global NUMERROR

[N,K]=size(A);
L=size(Y,1);
GradP=zeros(L,N);
R=sum(N - (1:(N-1)));

for k=1:K

    yk=Y(:,k);
    zk=Z(:,k);
    ak=A(:,k);
    byk=yk'*yk;
    dk=D(k);    
    Mk=diag((1-dk)*ones(L,1)+dk*zk);
    GradP=GradP-Mk'*yk*ak'/byk + Mk'*Mk*Po*ak*ak'/byk;
    
end
O = N*eye(N) - ones(N,N);  
GradP=GradP/K+rho*Po*O/R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Optimal Step in Update Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numG=rho*trace(GradP*O*Po'+Po*O*GradP')/R/2;
denG=rho*trace(GradP*O*GradP')/R;
for k=1:K
    yk=Y(:,k);
    zk=Z(:,k);
    ak=A(:,k);
    dk=D(k);  
    byk=yk'*yk;
    Mk=diag((1-dk)*ones(L,1)+dk*zk);
    T1=Mk*GradP*ak;
    numG=numG+T1'*Mk*(Po*ak-yk)/byk/K;
    denG=denG+T1'*T1/byk/K;
end
alpha=max([0, numG/denG]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Stepest Descent Update of End-members Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_est=Po-alpha*GradP;
P_est(P_est<0) = 0;
P_est(isnan(P_est))=0;
P_est(isinf(P_est))=0;
P=P_est./repmat(sum(P_est),L,1);

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
