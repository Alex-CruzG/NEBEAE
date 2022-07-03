%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Analysis of Cuprite Dataset of Nonlinear Unmixing --->> Table 5 and
% Figure 5
%
% DUCD
% June/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('CupriteS1_R188.mat');
load('groundTruth_Cuprite_nEnd12.mat');
Z=Y/max(Y(:));
n=12;
Po=M(slctBnds',1:n);
[L,K]=size(Z);
k=1:L;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Paremeters of NEBEAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initcond=8;             % Initial condition of end-members matrix
rho=0.1;               % Similarity weight in end-members estimation
lambda=0.2;             % Entropy weight for abundance estimation
epsilon=1e-3;
maxiter=50;
parallel=0;
normalization=1;
downsampling=0.25;       % Downsampling in end-members estimation
disp_iter=0;            % Display partial performance in BEAE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters of BSNHU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=0.5;
alpha=0.5;
maxG=1;
delta=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute EBEAE Methodologys 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Cuprite Dataset');
disp(['Number of end-members=' num2str(n)]);

%%
tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('NEBEAE');
paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,disp_iter];
[P1,A1,D1,G1,Zh1]=NEBEAE2(Z,n,paramvec);
Tnbeae=toc;
disp(['Estimation Error = ' num2str(norm(Zh1-Z,'fro')/norm(Z,'fro'))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P1))]);
disp(['End-members Error SAM =' num2str(errorSAM(Po,P1))]);
disp(['Computation time = ' num2str(Tnbeae)]);


%%
tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('Supervised MLM');
P_min=-100;
options=optimset('fmincon');
options = optimset(options,'Display','off','Algorithm','sqp','MaxFunEvals',50000,'TolFun',1e-10,'TolCon',1e-8,'GradObj','off');
Aa2=zeros(n+1,K); % The first p variables are the abundances, the p+1'th variable contains the P value
Zh2=zeros(L,K);
P2=VCA(Z,n);
Aini=pinv(P2)*Z;
Aini(Aini<0)=0;
Aini=Aini./repmat(sum(Aini),[n,1]);
for i=1:K
    a_ini=[Aini(:,i); 0.0]; % Initialize with linear unmixing results, P=0
    % Sum-to-one applies to abundances, not P. P is restricted to [P_min,1]
    Aa2(:,i) = fmincon(@(a) sum((Z(:,i)-model_MLM(a,P2)).^2), a_ini,[],[],[ones(1,n) 0],1,[zeros(n,1); P_min],ones(n+1,1),[],options);
    Zh2(:,i) = model_MLM(Aa2(:,i),P2);
end
Tsmlm=toc;
P2=P2./repmat(sum(P2),[L,1]);
A2=Aa2(1:n,:);
D2=Aa2(n+1,:);
disp(['Estimation Error = ' num2str(norm(Zh2-Z,'fro')/norm(Z,'fro'))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P2))]);
disp(['End-members Error SAM =' num2str(errorSAM(Po,P2))]);
disp(['Computation time = ' num2str(Tsmlm)]);

tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('Unsupervised MLM');
[P3,A3, D3,Zh3,~]=unmix(Z,n,maxiter);
Tumlm=toc;
P3=P3./repmat(sum(P3),[L,1]);
disp(['Estimation Error = ' num2str(norm(Zh3-Z,'fro')/norm(Z,'fro'))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P3))]);
disp(['End-members Error SAM =' num2str(errorSAM(Po,P2))]);
disp(['Computation time = ' num2str(Tumlm)]);

%%
tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('BSNHU');
Avca=vcaBSNHU(Z,'endmembers',n);
Svca= hyperFclsBSNHU(Z,Avca)';
[P4,B4,A4,~]=sparseBilinearFanUnmixing(Z',Avca,Svca,maxG,delta,q,[]);
G4=createB(A4);
Zh4=P4*A4'+ B4*G4';
Tbsnhu=toc;
P4=P4./repmat(sum(P4),[L,1]);
disp(['Estimation Error = ' num2str(norm(Zh4-Z,'fro')/norm(Z,'fro'))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P4))]);
disp(['End-members Error SAM =' num2str(errorSAM(Po,P2))]);
disp(['Computation time = ' num2str(Tbsnhu)]);


figure;
subplot(311)
hist(D1,50);
%xlim([0,1]);
title('Nonlinear Interaction Level by NEBEAE');
xlabel('d_k');
grid on;
subplot(312)
hist(D2,50);
%xlim([0,1]);
title('Nonlinear Interaction Level by Supervised-NSUMMM');
xlabel('d_k');
grid on;
subplot(313)
hist(D3,50);
%xlim([0,1]);
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
imagesc(reshape(D1,nRow,nCol), [min([D1;D2';D3]), max([D1;D2';D3])]); colormap;
title('(b) NEBEAE');
subplot(1,4,3);
imagesc(reshape(D2,nRow,nCol),[min([D1;D2';D3]), max([D1;D2';D3])]);  
title('(c) Supervised-MMMNSU');
subplot(1,4,4);
imagesc(reshape(D3,nRow,nCol), [min([D1;D2';D3]), max([D1;D2';D3])]);  
title('(d) Unsupervised-UNSUBMMM');
orient landscape


%%
AAnd=string(' & ');
EEnd=string(' \\');
ErrorCuprite=[strcat("$E_z$ & ",num2str(norm(Zh1-Z,'fro')/norm(Z,'fro')), AAnd, ...
    num2str(norm(Zh2-Z,'fro')/norm(Z,'fro')), AAnd, num2str(norm(Zh3-Z,'fro')/norm(Z,'fro')), EEnd);
    strcat("$E_p$ & ",num2str(errorendmembers(Po,P1)),  AAnd, num2str(errorendmembers(Po,P2)), AAnd, ...
    num2str(errorendmembers(Po,P3)), EEnd);
    strcat("$E_p$ & ",num2str(errorSAM(Po,P1)), AAnd, num2str(errorSAM(Po,P2)), AAnd, ...
    num2str(errorSAM(Po,P3)), EEnd);     
    strcat("Computational time & ",num2str(Tnbeae), AAnd, num2str(Tsmlm), AAnd, num2str(Tumlm), EEnd)];
disp(ErrorCuprite)