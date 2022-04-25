%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Unmixing with Nonlinear Models for Urban Dataset --> Figure 6 & 7
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

load('Urban_F210.mat');
load('end4_groundTruth.mat');
Ao=A;
Po=M;
Z=Y(SlectBands,:)./maxValue;
[L,K]=size(Z);
k=1:L;
n=size(A,1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Paremeters of EBEAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initcond=5;             % Initial condition of end-members matrix
rho=1;               % Similarity weight in end-members estimation
lambda=0.25;             % Entropy weight for abundance estimation
epsilon=1e-3;
maxiter=50;
parallel=0;
normalization=1;
downsampling=0.0;       % Downsampling in end-members estimation
disp_iter=0;            % Display partial performance in BEAE
%%
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
disp(['Abundances Error =' num2str(errorabundances(Ao,A))]);
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
P3=PPo;
for i=1:K
    a_ini=[A2(:,i); 0.0];
    Aa3(:,i) = fmincon(@(a) sum((Z(:,i)-model_MLM(a,P3)).^2), a_ini,[],[],[ones(1,n) 0],1,[zeros(n,1); P_min],ones(n+1,1),[],options);
    Zh3(:,i) = model_MLM(Aa3(:,i),P3);
end
Tsmlm=toc;
P3=P3./repmat(sum(P3),[L,1]);
A3=Aa3(1:n,:);
D3=Aa3(n+1,:);
disp(['Estimation Error = ' num2str(norm(Zh3-Z,'fro')/norm(Z,'fro'))]);
disp(['Abundances Error =' num2str(errorabundances(Ao,A3))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P3))]);
disp(['Computation time = ' num2str(Tsmlm)]);

tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('Unsupervised MLM');
[P4,A4,D4,Zh4,~]=unmix(Z,n,maxiter);
Tumlm=toc;
P4=P4./repmat(sum(P4),[L,1]);
disp(['Estimation Error = ' num2str(norm(Zh4-Z,'fro')/norm(Z,'fro'))]);
disp(['Abundances Error =' num2str(errorabundances(Ao,A4))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P4))]);
disp(['Computation time = ' num2str(Tumlm)]);


%%
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



figure;
subplot(411);
imagesc([reshape(Ao(1,:),nRow,nCol) reshape(Ao(2,:),nRow,nCol)... 
       reshape(Ao(3,:),nRow,nCol) reshape(Ao(4,:),nRow,nCol)]);
title('Ground-truth abundances');
subplot(412);
%A1=A./sum(A);
imagesc([reshape(A1(1,:),nRow,nCol) reshape(A1(2,:),nRow,nCol)... 
    reshape(A1(3,:),nRow,nCol) reshape(A1(4,:),nRow,nCol)]);
title('Estimated abundances by NEBEAE');
subplot(413);
imagesc([reshape(A3(1,:),nRow,nCol) reshape(A3(2,:),nRow,nCol)... 
    reshape(A3(3,:),nRow,nCol) reshape(A3(4,:),nRow,nCol)]);
title('Estimated abundances by Supervised MLM');
subplot(414);
imagesc([reshape(A4(1,:),nRow,nCol) reshape(A4(2,:),nRow,nCol)... 
    reshape(A4(3,:),nRow,nCol) reshape(A4(4,:),nRow,nCol)]);
title('Estimated abundances by Unsupervised MLM');

RGBU=imread('UrbanRGB.png');

figure;
subplot(1,5,1);
imagesc(RGBU);
axis off
title({'(a) Aproximated','RGB Image'});
subplot(1,5,2);
imagesc(reshape(Ao(1,:),nRow,nCol))
title({'(b) Abundance','End-member 1'});
subplot(1,5,3);
imagesc(reshape(Ao(2,:),nRow,nCol));
title({'(c) Abundance','End-member 2'});
subplot(1,5,4);
imagesc(reshape(Ao(3,:),nRow,nCol));
title({'(d) Abundance','End-member 3'});
subplot(1,5,5);
imagesc(reshape(Ao(4,:),nRow,nCol));
title({'(e) Abundance','End-member 4'});
orient landscape

%%
figure
imagesc(reshape(D,nRow,nCol))
title('Nonlinear Interaction Level');
colorbar


