%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Synthetic Evaluation of Unmixing with Nonlinear Model --> Figure 4
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

SNR=20;     
PSNR=20;
Nsamples=60;
n=4;
ModelType=5;

[Z,Po,Ao,Go]=VNIRsynth(n,Nsamples,SNR,PSNR,ModelType);
[L,K]=size(Z);
Y=Z;
k=1:L;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Paremeters of EBEAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initcond=6;             % Initial condition of end-members matrix
rho=0.1;               % Similarity weight in end-members estimation
lambda=0.1;             % Entropy weight for abundance estimation
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
disp(['SNR =' num2str(SNR) ' dB']);
disp(['PSNR =' num2str(PSNR) ' dB']);
disp(['Number of end-members=' num2str(n)]);

%%
tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('NEBEAE');
paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,disp_iter];
[P,A,D,G,Zh]=NEBEAE(Z,n,paramvec);
Tnbeae1=toc;
disp(['Estimation Error = ' num2str(norm(Zh-Z,'fro')/norm(Z,'fro'))]);
disp(['Abundances Error =' num2str(errorabundances(Ao,(A.*repmat(G',[n,1]))))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P))]);
disp(['Computation time = ' num2str(Tnbeae1)]);

%%
%%
tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('Supervised MLM');
P_min=-100;
options=optimset('fmincon');
options = optimset(options,'Display','off','Algorithm','sqp','MaxFunEvals',50000,'TolFun',1e-10,'TolCon',1e-8,'GradObj','off');
Aa3=zeros(n+1,K); % The first p variables are the abundances, the p+1'th variable contains the P value
Zh3=zeros(L,K);
P3=VCA(Z,n);
A2=pinv(P3)*Z;
for i=1:K
    a_ini=[A2(:,i); 0.0]; % Initialize with linear unmixing results, P=0
    % Sum-to-one applies to abundances, not P. P is restricted to [P_min,1]
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
[P4,A4, D4,Zh4,~]=unmix(Z,n,maxiter);
Tumlm=toc;
P4=P4./repmat(sum(P4),[L,1]);
disp(['Estimation Error = ' num2str(norm(Zh4-Z,'fro')/norm(Z,'fro'))]);
disp(['Abundances Error =' num2str(errorabundances(Ao,A4))]);
disp(['End-members Error =' num2str(errorendmembers(Po,P4))]);
disp(['Computation time = ' num2str(Tumlm)]);



%%


figure;
subplot(3,5,5)
hist(D);
ylabel('histogram');
xlim([0,1]);
title({'Nonlinear Interaction','Level'});
xlabel('d_k');
grid on;
subplot(3,5,10)
hist(D3);
xlim([0,1]);
ylabel('histogram');
title({'Nonlinear Interaction','Level'});
xlabel('d_k');
grid on;
subplot(3,5,15)
hist(D4);
ylabel('histogram');
xlim([0,1]);
title({'Nonlinear Interaction','Level'});
xlabel('d_k');
grid on;

subplot(3,5,1);
imagesc(reshape(A(1,:),Nsamples,Nsamples))
title({'Abundance','End-member 1'})
subplot(3,5,2);
imagesc(reshape(A(2,:),Nsamples,Nsamples))
title({'Abundance','End-member 2'})
subplot(3,5,3);
imagesc(reshape(A(3,:),Nsamples,Nsamples))
title({'Abundance','End-member 3'})
subplot(3,5,4);
imagesc(reshape(A(4,:),Nsamples,Nsamples));
title({'Abundance','End-member 4'})
subplot(3,5,6);
imagesc(reshape(A3(1,:),Nsamples,Nsamples))
title({'Abundance','End-member 1'})
subplot(3,5,7);
imagesc(reshape(A3(2,:),Nsamples,Nsamples))
title({'Abundance','End-member 2'})
subplot(3,5,8);
imagesc(reshape(A3(3,:),Nsamples,Nsamples))
title({'Abundance','End-member 3'})
subplot(3,5,9);
imagesc(reshape(A3(4,:),Nsamples,Nsamples));
title({'Abundance','End-member 4'})
subplot(3,5,11);
imagesc(reshape(A4(1,:),Nsamples,Nsamples))
title({'Abundance','End-member 1'})
subplot(3,5,12);
imagesc(reshape(A4(2,:),Nsamples,Nsamples))
title({'Abundance','End-member 2'})
subplot(3,5,13);
imagesc(reshape(A4(3,:),Nsamples,Nsamples)) 
title({'Abundance','End-member 3'})
subplot(3,5,14);
imagesc(reshape(A4(4,:),Nsamples,Nsamples));
title({'Abundance','End-member 4'})
orient landscape





