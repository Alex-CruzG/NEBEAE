%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Unmixing with Nonlinear Models for Curprite Dataset --> Figure 5
%
% ``Nonlinear Extended Blind End-member and Abundance Extraction for
% Hyperspectral Images''
%  Campos-Delgado D.U. et al, Submitted to Signal Processing (Elsevier)
%
%
% DUCD
% April/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('groundTruth_Cuprite_nEnd12.mat');
load('CupriteS1_R188.mat');
Z=Y/max(Y(:));
n=12;
Po=M(slctBnds',1:n);
[L,K]=size(Z);
k=1:L;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Paremeters of EBEAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initcond=5;             % Initial condition of end-members matrix
rho=0.1;               % Similarity weight in end-members estimation
lambda=0.2;             % Entropy weight for abundance estimation
epsilon=1e-3;
maxiter=50;
parallel=0;
normalization=1;
downsampling=0.25;       % Downsampling in end-members estimation
disp_iter=0;            % Display partial performance in BEAE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Nonlinear Unmixing Methodologies 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Synthetic VNIR Datasets');
disp('Blind Unmixing Estimation');
disp(['Number of end-members=' num2str(n)]);

tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('NEBEAE');
paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,disp_iter];
[P,A,D,G,Zh]=NEBEAE(Z,n,paramvec);
Tnbeae1=toc;
disp(['Estimation Error = ' num2str(norm(Zh-Z,'fro')/norm(Z,'fro'))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P))]);
disp(['Computation time = ' num2str(Tnbeae1)]);


%%
tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('Supervised MLM');
P_min=-100;
options=optimset('fmincon');
options = optimset(options,'Display','off','Algorithm','sqp','MaxFunEvals',50000,'TolFun',1e-10,'TolCon',1e-8,'GradObj','off');
Aa3=zeros(n+1,K); 
Zh3=zeros(L,K);
P3=VCA(Z,n);
A1=pinv(P3)*Z;
for i=1:K
    a_ini=[A1(:,i); 0.0]; 
    Aa3(:,i) = fmincon(@(a) sum((Z(:,i)-model_MLM(a,P3)).^2), a_ini,[],[],[ones(1,n) 0],1,[zeros(n,1); P_min],ones(n+1,1),[],options);
    Zh3(:,i) = model_MLM(Aa3(:,i),P3);
end
Tsmlm=toc;
P3=P3./repmat(sum(P3),[L,1]);
A3=Aa3(1:n,:);
D3=Aa3(n+1,:);
disp(['Estimation Error = ' num2str(norm(Zh3-Z,'fro')/norm(Z,'fro'))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P3))]);
disp(['Computation time = ' num2str(Tsmlm)]);

tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('Unsupervised MLM');
[P4,A4, D4,Zh4,~]=unmix(Z,n,maxiter);
Tumlm=toc;
P4=P4./repmat(sum(P4),[L,1]);
disp(['Estimation Error = ' num2str(norm(Zh4-Z,'fro')/norm(Z,'fro'))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P4))]);
disp(['Computation time = ' num2str(Tumlm)]);

figure;
subplot(311)
hist(D,50);
title('Nonlinear Interaction Level by NEBEAE');
xlabel('d_k');
grid on;
subplot(312)
hist(D3,50);
title('Nonlinear Interaction Level by Supervised-NSUMMM');
xlabel('d_k');
grid on;
subplot(313)
hist(D4,50);
title('Nonlinear Interaction Level by Unsupervised-NSUMMM');
xlabel('d_k');
grid on;
%%

CRGB=imread('CupriteRGB2.png');

figure;
subplot(1,4,1);
imagesc(CRGB);
axis off
title('(a) Approximated RGB Image');
subplot(1,4,2);
imagesc(reshape(D,nRow,nCol), [min([D;D3';D4]), max([D;D3';D4])]); colormap;
title('(b) NEBEAE');
subplot(1,4,3);
imagesc(reshape(D3,nRow,nCol),[min([D;D3';D4]), max([D;D3';D4])]);  
title('(c) Supervised-NSUMMM');
subplot(1,4,4);
imagesc(reshape(D4,nRow,nCol), [min([D;D3';D4]), max([D;D3';D4])]);  
title('(d) Unsupervised-NSUMMM');
orient landscape

%%
ErrorCuprita=[norm(Zh-Z,'fro')/norm(Z,'fro') norm(Zh3-Z,'fro')/norm(Z,'fro') norm(Zh4-Z,'fro')/norm(Z,'fro');
        errorendmembers(Po,P)  errorendmembers(Po,P3) errorendmembers(Po,P4);
        Tnbeae1 Tsmlm Tumlm];
