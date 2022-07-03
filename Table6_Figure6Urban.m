%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Analysis of Urban Dataset of Nonlinear Unmixing --->> Table 6 and
% Figure 6
%
% DUCD
% June/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Urban_F210.mat');
ErrorUrbanZ=zeros(3,4);
ErrorUrbanA=zeros(3,4);
ErrorUrbanP=zeros(3,4);
ErrorUrbanSAM=zeros(3,4);
ErrorUrbanT=zeros(3,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Paremeters of NEBEAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initcond=5;             % Initial condition of end-members matrix
rho=1;               % Similarity weight in end-members estimation
lambda=0.25;             % Entropy weight for abundance estimation
epsilon=1e-3;
maxiter=50;
parallel=0;
normalization=1;
downsampling=0.0;       % Downsampling in end-members estimation
disp_iter=0;            % Display partial performance in BEAE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters of BSNHU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q=0.5;
alpha=0.5;
maxG=1;
delta=1;

for kk=1:3
    if kk==1
        load('end4_groundTruth.mat');
    elseif kk==2
        load('end5_groundTruth.mat');
    else
        load('end6_groundTruth.mat');
    end
    Ao=A;
    Po=M;
    Z=Y(SlectBands,:)./maxValue;
    [L,K]=size(Z);
    k=1:L;
    n=size(A,1);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Execute NEBEAE Methodology
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('Urban Dataset');
    disp('Blind Unmixing Estimation');
    disp(['Number of end-members=' num2str(n)]);  

%%
    tic;
    disp('%%%%%%%%%%%%%%%%%%');
    disp('NEBEAE');
    paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,disp_iter];
    [P1,A1,D1,G1,Zh1]=NEBEAE2(Z,n,paramvec);
    Tnbeae=toc;
    ErrorUrbanZ(kk,1)=norm(Zh1-Z,'fro')/norm(Z,'fro');
    ErrorUrbanA(kk,1)=errorabundances(Ao,A1);
    ErrorUrbanP(kk,1)=errorendmembers(Po,P1);
    ErrorUrbanSAM(kk,1)=errorSAM(Po,P1);
    ErrorUrbanT(kk,1)=Tnbeae;
    disp(['Estimation Error = ' num2str(ErrorUrbanZ(kk,1))]);
    disp(['Abundances Error ='  num2str(ErrorUrbanA(kk,1))]);
    disp(['End-members Error =' num2str(ErrorUrbanP(kk,1))]);
    disp(['End-members SAM =' num2str(ErrorUrbanSAM(kk,1))]);
    disp(['Computation time = ' num2str(ErrorUrbanT(kk,1))]);

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
    Aini=pinv(P3)*Z;
    Aini(Aini<0)=0;
    Aini=Aini./repmat(sum(Aini),[n,1]);
    for i=1:K
        a_ini=[Aini(:,i); 0.0]; % Initialize with linear unmixing results, P=0
        % Sum-to-one applies to abundances, not P. P is restricted to [P_min,1]
        Aa3(:,i) = fmincon(@(a) sum((Z(:,i)-model_MLM(a,P3)).^2), a_ini,[],[],[ones(1,n) 0],1,[zeros(n,1); P_min],ones(n+1,1),[],options);
        Zh3(:,i) = model_MLM(Aa3(:,i),P3);
    end
    Tsmlm=toc;
    P3=P3./repmat(sum(P3),[L,1]);
    A3=Aa3(1:n,:);
    D3=Aa3(n+1,:);
    ErrorUrbanZ(kk,2)=norm(Zh3-Z,'fro')/norm(Z,'fro');
    ErrorUrbanA(kk,2)=errorabundances(Ao,A3);
    ErrorUrbanP(kk,2)=errorendmembers(Po,P3);
    ErrorUrbanSAM(kk,2)=errorSAM(Po,P3);
    ErrorUrbanT(kk,2)=Tsmlm;
    disp(['Estimation Error = ' num2str(ErrorUrbanZ(kk,2))]);
    disp(['Abundances Error ='  num2str(ErrorUrbanA(kk,2))]);
    disp(['End-members Error =' num2str(ErrorUrbanP(kk,2))]);
    disp(['End-members SAM =' num2str(ErrorUrbanSAM(kk,2))]);
    disp(['Computation time = ' num2str(ErrorUrbanT(kk,2))]);

    tic;
    disp('%%%%%%%%%%%%%%%%%%');
    disp('Unsupervised MLM');
    [P4,A4,D4,Zh4,~]=unmix(Z,n,maxiter);
    Tumlm=toc;
    P4=P4./repmat(sum(P4),[L,1]);
    ErrorUrbanZ(kk,3)=norm(Zh4-Z,'fro')/norm(Z,'fro');
    ErrorUrbanA(kk,3)=errorabundances(Ao,A4);
    ErrorUrbanP(kk,3)=errorendmembers(Po,P4);
    ErrorUrbanSAM(kk,3)=errorSAM(Po,P4);
    ErrorUrbanT(kk,3)=Tumlm;
    disp(['Estimation Error = ' num2str(ErrorUrbanZ(kk,3))]);
    disp(['Abundances Error ='  num2str(ErrorUrbanA(kk,3))]);
    disp(['End-members Error =' num2str(ErrorUrbanP(kk,3))]);
    disp(['End-members SAM =' num2str(ErrorUrbanSAM(kk,3))]);
    disp(['Computation time = ' num2str(ErrorUrbanT(kk,3))]);
    
    tic;
    disp('%%%%%%%%%%%%%%%%%%');
    disp('BSNHU');
    Avca=vcaBSNHU(Z,'endmembers',n);
    Svca= hyperFclsBSNHU(Z,Avca)';
    [P5,B5,A5,~]=sparseBilinearFanUnmixing(Z',Avca,Svca,maxG,delta,q,[]);
    G5=createB(A5);
    Zh5=P5*A5'+ B5*G5';
    Tbsnhu=toc;
    P5=P5./repmat(sum(P5),[L,1]);
    ErrorUrbanZ(kk,4)=norm(Zh5-Z,'fro')/norm(Z,'fro');
    ErrorUrbanA(kk,4)=errorabundances(Ao,A5);
    ErrorUrbanP(kk,4)=errorendmembers(Po,P5);
    ErrorUrbanSAM(kk,4)=errorSAM(Po,P5);
    ErrorUrbanT(kk,4)=Tbsnhu;
    disp(['Estimation Error = ' num2str(ErrorUrbanZ(kk,4))]);
    disp(['Abundances Error ='  num2str(ErrorUrbanA(kk,4))]);
    disp(['End-members Error =' num2str(ErrorUrbanP(kk,4))]);
    disp(['End-members SAM =' num2str(ErrorUrbanSAM(kk,4))]);
    disp(['Computation time = ' num2str(ErrorUrbanT(kk,4))]);

 end
  %%
 AAnd=string(' & ');
 EEnd=string(' \\');
 display('Estimation Error');
 ErrorUz=[strcat("4 & ",num2str(ErrorUrbanZ(1,1)), AAnd, ...
    num2str(ErrorUrbanZ(1,2)), AAnd, num2str(ErrorUrbanZ(1,3)), AAnd, num2str(ErrorUrbanZ(1,4)), EEnd);
    strcat("5 & ",num2str(ErrorUrbanZ(2,1)), AAnd, ...
    num2str(ErrorUrbanZ(2,2)), AAnd, num2str(ErrorUrbanZ(2,3)), AAnd, num2str(ErrorUrbanZ(2,4)), EEnd);
    strcat("6 & ",num2str(ErrorUrbanZ(3,1)), AAnd, ...
    num2str(ErrorUrbanZ(3,2)), AAnd, num2str(ErrorUrbanZ(3,3)), AAnd, num2str(ErrorUrbanZ(3,4)), EEnd)];
 display(ErrorUz);
 display('Abundance Error');
  ErrorUa=[strcat("4 & ",num2str(ErrorUrbanA(1,1)), AAnd, ...
    num2str(ErrorUrbanA(1,2)), AAnd, num2str(ErrorUrbanA(1,3)), AAnd, num2str(ErrorUrbanA(1,4)), EEnd);
    strcat("5 & ",num2str(ErrorUrbanA(2,1)), AAnd, ...
    num2str(ErrorUrbanA(2,2)), AAnd, num2str(ErrorUrbanA(2,3)), AAnd, num2str(ErrorUrbanA(2,4)), EEnd);
    strcat("6 & ",num2str(ErrorUrbanA(3,1)), AAnd, ...
    num2str(ErrorUrbanA(3,2)), AAnd, num2str(ErrorUrbanA(3,3)), AAnd, num2str(ErrorUrbanA(3,4)), EEnd)];
 display(ErrorUa);
 display('End-member Error');
   ErrorUp=[strcat("4 & ",num2str(ErrorUrbanP(1,1)), AAnd, ...
    num2str(ErrorUrbanP(1,2)), AAnd, num2str(ErrorUrbanP(1,3)), AAnd, num2str(ErrorUrbanP(1,4)), EEnd);
    strcat("5 & ",num2str(ErrorUrbanP(2,1)), AAnd, ...
    num2str(ErrorUrbanP(2,2)), AAnd, num2str(ErrorUrbanP(2,3)), AAnd, num2str(ErrorUrbanP(2,4)), EEnd);
    strcat("6 & ",num2str(ErrorUrbanP(3,1)), AAnd, ...
    num2str(ErrorUrbanP(3,2)), AAnd, num2str(ErrorUrbanP(3,3)), AAnd, num2str(ErrorUrbanP(3,4)), EEnd)];
 display(ErrorUp);
 display('End-member Error SAM');
    ErrorUpsam=[strcat("4 & ",num2str(ErrorUrbanSAM(1,1)), AAnd, ...
    num2str(ErrorUrbanSAM(1,2)), AAnd, num2str(ErrorUrbanSAM(1,3)), AAnd, num2str(ErrorUrbanSAM(1,4)), EEnd);
    strcat("5 & ",num2str(ErrorUrbanSAM(2,1)), AAnd, ...
    num2str(ErrorUrbanSAM(2,2)), AAnd, num2str(ErrorUrbanSAM(2,3)), AAnd, num2str(ErrorUrbanSAM(2,4)), EEnd);
    strcat("6 & ",num2str(ErrorUrbanSAM(3,1)), AAnd, ...
    num2str(ErrorUrbanSAM(3,2)), AAnd, num2str(ErrorUrbanSAM(3,3)), AAnd, num2str(ErrorUrbanSAM(3,4)), EEnd)];
 display(ErrorUpsam);
 display('Computation Time');
     ErrorUt=[strcat("4 & ",num2str(ErrorUrbanT(1,1)), AAnd, ...
    num2str(ErrorUrbanT(1,2)), AAnd, num2str(ErrorUrbanT(1,3)), AAnd, num2str(ErrorUrbanT(1,4)), EEnd);
    strcat("5 & ",num2str(ErrorUrbanT(2,1)), AAnd, ...
    num2str(ErrorUrbanT(2,2)), AAnd, num2str(ErrorUrbanT(2,3)), AAnd, num2str(ErrorUrbanT(2,4)), EEnd);
    strcat("6 & ",num2str(ErrorUrbanT(3,1)), AAnd, ...
    num2str(ErrorUrbanT(3,2)), AAnd, num2str(ErrorUrbanT(3,3)), AAnd, num2str(ErrorUrbanT(3,4)), EEnd)];
 display(ErrorUt);
  
    

%%
figure;
subplot(311)
hist(D1,50);
%xlim([0,1]);
title('Nonlinear Interaction Level by NEBEAE');
xlabel('d_k');
grid on;
subplot(312)
hist(D3,50);
%xlim([0,1]);
title('Nonlinear Interaction Level by Supervised-MMMNSU');
xlabel('d_k');
grid on;
subplot(313)
hist(D4,50);
%xlim([0,1]);
title('Nonlinear Interaction Level by Unsupervised-UNSUBMMM');
xlabel('d_k');
grid on;
%%
UrbanRGB=imread('UrbanRGB.png');

figure;
subplot(1,4,1);
imagesc(UrbanRGB);
axis off
title('(a) Approximated RGB Image');
subplot(1,4,2);
imagesc(reshape(D1,nRow,nCol), [min([D1;D3';D4]), max([D1;D3';D4])]); colormap;
title('(b) NEBEAE');
subplot(1,4,3);
imagesc(reshape(D3,nRow,nCol),[min([D1;D3';D4]), max([D1;D3';D4])]);  
title('(c) Supervised-MMMNSU');
subplot(1,4,4);
imagesc(reshape(D4,nRow,nCol), [min([D1;D3';D4]), max([D1;D3';D4])]);  
title('(d) Unsupervised-UNSUBMMM');
orient landscape

