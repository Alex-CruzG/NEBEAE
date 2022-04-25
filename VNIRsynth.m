function [Y,P,A,G]=VNIRsynth(N,Npixels,SNR,PSNR,ModelType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [Y,P,A,G,Yo,D]=VNIRsynth2(N,Npixels,SNR,PSNR,ModelType)
%
% Synthetic mFLIM dataset with linear or non-linear mixing models:
%
% 0) Linear Mixing Model (LMM)
%
% y_k = \sum_{n=1}^N a_{k,n} p_n + v_k         v_k \sim N(0,\sigma^2 I)
%
% 1) Fan Model (FM) --> Fan et al., 2009
%
% y_k = \sum_{n=1}^N a_{k,n} p_n 
%         + \sum_{n=1}^{N-1} \sum_{m=n+1}^N (a_{k,n}p_n) \odot (a_{k,m}p_m) + v_k
%
% 2) Generalized Bilinear Model (GBM) --> Halimi et al., 2011)
%
% y_k = \sum_{n=1}^N a_{k,n} p_n 
%         + \sum_{n=1}^{N-1} \sum_{m=1}^N \gamma_{n,m} (a_{k,n}p_n) \odot (a_{k,m}p_m) + v_k
%
% 3) Postnonlinear Mixing Model (PNMM) --> Altmann et al., 2012
%
% y_k = \sum_{n=1}^N a_{k,n} p_n 
%         + \sum_{n=1}^{N} \sum_{m=1}^N \xi (a_{k,n}p_n) \odot (a_{k,m}p_m) + v_k
%
% 4) Multilinear Mixing Model (MMM) --> Heylen and Scheunders (2016)
%
%  y_k = (1-P_k) \sum_{n=1}^N a_{k,n} p_n / (1-P_k \sum_{n=1}^N a_{k,n} p_n) + v_k
%
% INPUTS
% N --> Order of multi-exponential model \in [2,4]
% Npixels --> numbers of pixels in x & y axes
% SNR --> SNR of Gaussian noise (dB)
% PSNR --> PSNR for Shot noise (dB)
% ModelType --> 0 (LMM-Default), 1 (FM), 2 (GBM), 3 (PNMM) and 4 (MMM)
%
% OUTPUTS
% Y --> matrix of fluorescence decays of size 186 x (Npixels*Npixels)
% P --> matrix of end-members 186 x N
% A --> matrix of abundances of N x (Npixels*Npixels)
% G --> matrix of nonlinear interaction 186 x (Npixels*Npixels)
% Yo --> noiseless measurements
% Do --> nonlinear interaction level
%
% Y = P*A+G
%
% Feb/2021
% DUCD
%
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthetic VNIR Dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    N=2;
    Npixels=60;
    SNR=40;
    PSNR=20;
    ModelType=0;
elseif nargin<5
    ModelType=0;
elseif nargin<4
    PSNR=0;
    ModelType=0;
elseif nargin<3
    PSNR=0; SNR=0; ModelType=0;
elseif nargin<2 
    PSNR=0; SNR=0; Npixels=60;ModelType=0;
end

if N>4,
    N=4; disp('The maximum number of components is 4!!');
end

if SNR ~= 0,
   NoiseGaussian=1;
else
   NoiseGaussian=0;
end
if PSNR ~= 0,
   NoiseShot=1; 
else
    NoiseShot=0;
end
if SNR ~= 0 || PSNR ~= 0,
    NoiseMeasurement=1;        
else
     NoiseMeasurement=0;
end

Nsamp=Npixels;
x=1:Npixels;
y=1:Npixels;
[xx,yy]=meshgrid(x,y);
K=Npixels*Npixels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End-members extracted from
% https://crustal.usgs.gov/speclab/QueryAll07a.php
% Kokaly, R.F. et al., 2017, 
% USGS Spectral Library Version 7: U.S. Geological Survey Data Series 1035, 61 p., 
% https://doi.org/10.3133/ds1035.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%AnthophylliteWLf;AnthophylliteEMf;  
%AndalusiteWLf; AndalusiteEMf;
BloediteWLf; BloediteEMf;
BluteriteEMf; BluteriteWLf;
AxiniteWLf;AxiniteEMf;
%AlmandineWLf; AlmandineEMf;
%AlbiteWLf; AlbiteEMf;
AntigoriteEMf; AntigoriteWLf;
CarbonEMf; CarbonWLf;

%WL1=CarbonWL(CarbonEM>0); EM1=CarbonEM(CarbonEM>0);
WL2=BloediteWL(BloediteEM>0); EM2=BloediteEM(BloediteEM>0); 
WL3=BluteriteWL(BluteriteEM>0); EM3=BluteriteEM(BluteriteEM>0); 
WL4=AntigoriteWL(AntigoriteEM>0); EM4=AntigoriteEM(AntigoriteEM>0);

%WL1=AnthophylliteWL(AnthophylliteEM>0); EM1=AnthophylliteEM(AnthophylliteEM>0);
%WL2=AndalusiteWL(AndalusiteEM>0); EM2=AndalusiteEM(AndalusiteEM>0);
WL1=AxiniteWL(AxiniteEM>0); EM1=AxiniteEM(AxiniteEM>0);
%WL3=AlmandineWL(AlmandineEM>0); EM3=AlmandineEM(AlmandineEM>0);
%WL4=AlbiteWL(AlbiteEM>0); EM4=AlbiteEM(AlbiteEM>0);

L=size(WL4,1);
%WL=WL1;
EM=[EM1(1:L) EM2(1:L) EM3(1:L) EM4(1:L)];
P1=EM1(1:L)/max(max(EM));
P2=EM2(1:L)/max(max(EM));
P3=EM3(1:L)/max(max(EM));
P4=EM4(1:L)/max(max(EM));

gamma0=[0.5 0.3 0.25 0.5 0.6 0.2]; % (2) Mixing coefficients in GBM
xi0=0.3;                          % (3) Scaling coefficient in PNMM
prob=ones(Npixels,Npixels)*0.5;                          % (4) Probability of nonlinear mixing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Abundance Maps and Nonlinear Interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if N==2

    aa1=7*exp(-0.001*(xx-Nsamp/2).^2-0.001*(yy-Nsamp/2).^2)+1.0;
    aa2=2.5*exp(-0.001*(xx-Nsamp).^2-0.001*(yy).^2)+2.5*exp(-0.0001*(xx).^2-0.001*(yy-Nsamp).^2)+...
            2.5*exp(-0.001*(xx).^2-0.0001*(yy-Nsamp).^2)+2.5*exp(-0.0001*(xx-Nsamp).^2-0.001*(yy-Nsamp).^2);

     a1=zeros(Nsamp,Nsamp);
     a2=zeros(Nsamp,Nsamp);
    for i=1:Nsamp,
        for l=1:Nsamp,
            a1(i,l)=aa1(i,l)/(aa1(i,l)+aa2(i,l));
            a2(i,l)=aa2(i,l)/(aa1(i,l)+aa2(i,l));
        end
    end

elseif N==3,

    aa1=7*exp(-0.005*(xx-Nsamp/2).^2-0.005*(yy-Nsamp/2).^2)+0.5;
    aa2=2.5*exp(-0.001*(xx-Nsamp).^2-0.001*(yy).^2)+2.5*exp(-0.0001*(xx).^2-0.001*(yy-Nsamp).^2);
    aa3=3.5*exp(-0.001*(xx).^2-0.0001*(yy-Nsamp).^2)+2.5*exp(-0.0001*(xx-Nsamp).^2-0.001*(yy-Nsamp).^2);

    a1=zeros(Nsamp,Nsamp);
    a2=zeros(Nsamp,Nsamp);
    a3=zeros(Nsamp,Nsamp);

    for i=1:Nsamp,
        for l=1:Nsamp,
            a1(i,l)=aa1(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l));
            a2(i,l)=aa2(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l));
            a3(i,l)=aa3(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l));
        end
    end

  
elseif N==4,
    
    aa1=2.5*exp(-0.005*(xx-Nsamp/2).^2-0.0005*(yy-Nsamp/2).^2)+0.5;
    aa2=2.5*exp(-0.001*(xx-Nsamp).^2-0.00025*(yy).^2);%+2.5*exp(-0.0001*(xx).^2-0.001*(yy-Nsamp).^2);
    aa3=2.5*exp(-0.001*(xx).^2-0.0002*(yy-Nsamp).^2);
    aa4=2.5*exp(-0.001*(xx-8*Nsamp/9).^2-0.001*(yy-8*Nsamp/9).^2)+2.5*exp(-0.001*(xx-Nsamp/9).^2-0.001*(yy-8*Nsamp/9).^2);
    
    a1=zeros(Nsamp,Nsamp);
    a2=zeros(Nsamp,Nsamp);
    a3=zeros(Nsamp,Nsamp);
    a4=zeros(Nsamp,Nsamp);
    for i=1:Nsamp
        for l=1:Nsamp
            a1(i,l)=aa1(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a2(i,l)=aa2(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a3(i,l)=aa3(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a4(i,l)=aa4(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
        end
    end

end

Yy=zeros(Nsamp,Nsamp,L);
Gg=zeros(Nsamp,Nsamp,L);
Dd=zeros(Nsamp,Nsamp);

for i=1:Nsamp
    for j=1:Nsamp
        if N==2
            if ModelType==1
                g=(a1(i,j)*P1).*(a2(i,j)*P2);
            elseif ModelType==2
                %gamma=rand(sum(1:(N-1)),1);
                g=(a1(i,j)*P1).*(a2(i,j)*P2)*gamma(1);
            elseif ModelType==3
                g=(a1(i,j)*P1 + a2(i,j)*P2).*(a1(i,j)*P1 + a2(i,j)*P2)*xi;
            else
                g=0;
            end
            y=a1(i,j)*P1 + a2(i,j)*P2;
        elseif N==3
            if ModelType==1
                 g=(a1(i,j)*P1).*(a2(i,j)*P2) + (a1(i,j)*P1).*(a3(i,j)*P3) + (a2(i,j)*P2).*(a3(i,j)*P3);
            elseif ModelType==2
                %gamma=rand(sum(1:(N-1)),1);
                gammaL=gamma0+randn(6)*0.1;
                 g=(a1(i,j)*P1).*(a2(i,j)*P2)*gammaL(1)+(a1(i,j)*P1).*(a3(i,j)*P3)*gammaL(2)+(a2(i,j)*P2).*(a3(i,j)*P3)*gammaL(3);
            elseif ModelType==3
                 xi=xi0+randn*0.01;
                 g=(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3).*(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3)*xi;
            else
                g=0;
            end
            y=a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3; 
        elseif N==4
            if ModelType==1
                g1=(a1(i,j)*P1).*(a2(i,j)*P2) + (a1(i,j)*P1).*(a3(i,j)*P3) + (a1(i,j)*P1).*(a4(i,j)*P4);
                g=g1 + (a2(i,j)*P2).*(a3(i,j)*P3) + (a2(i,j)*P2).*(a4(i,j)*P4) + (a3(i,j)*P3).*(a4(i,j)*P4);
            elseif ModelType==2
                %gamma=rand(sum(1:(N-1)),1);
                g1=(a1(i,j)*P1).*(a2(i,j)*P2)*gamma(1) + (a1(i,j)*P1).*(a3(i,j)*P3)*gamma(2) + (a1(i,j)*P1).*(a4(i,j)*P4)*gamma(3);
                g=g1 + (a2(i,j)*P2).*(a3(i,j)*P3)*gamma(4) + (a2(i,j)*P2).*(a4(i,j)*P4)*gamma(5) + (a3(i,j)*P3).*(a4(i,j)*P4)*gamma(6);
            elseif ModelType==3
                g=(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3 + a4(i,j)*P4).*(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3 + a4(i,j)*P4)*xi;
            else
                g=0;
            end
            y=a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3 + a4(i,j)*P4; 
        end
                      
        Gg(i,j,:)=g;


        if ModelType==5
            %d=max([0,abs(Pk*randn)]);
            x=y./sum(y);
            y=((1-prob(i,j))*x)./(1-(prob(i,j)*x));
            Dd(i,j)=prob(i,j);
        else 
            Dd(i,j)=0;
        end
        Gg(i,j,:)=g;     
        Yy(i,j,:)=y+g;
   
    end
end

Yo=reshape(Yy,K,L)';
Do=reshape(Dd,K,1);
%Ym=max(reshape(Yy,K,L)',[],2);
Ym=mean(reshape(Yy,K,L)',2);

if NoiseMeasurement==1 && NoiseGaussian==1
    sigmay=sqrt((1/(L-1))*(Ym'*Ym)/(10^(SNR/10)));
    Yy=Yy+sigmay*randn(Nsamp,Nsamp,L);
end

if NoiseMeasurement==1 && NoiseShot==1
     sigmay=sqrt(max(Ym)^2 /(10^(PSNR/10)));
     y1=(sigmay*repmat(Ym,1,Nsamp*Nsamp).*randn(L,Nsamp*Nsamp))';
     Yy=Yy+reshape(y1,Nsamp,Nsamp,L);
end  

if N==2     
    P=[P1 P2];
    A=[reshape(a1,1,K);reshape(a2,1,K)];
elseif N==3
    P=[P1 P2 P3];
    A=[reshape(a1,1,K);reshape(a2,1,K); reshape(a3,1,K)];
elseif N==4
    P=[P1 P2 P3 P4];
    A=[reshape(a1,1,K);reshape(a2,1,K); reshape(a3,1,K); reshape(a4,1,K)];
end
G=reshape(Gg,K,L)';
%Pm=sum(P);
%P=P./repmat(Pm,[L,1]);
Y=reshape(Yy,K,L)';
%Y=Y/max(Y(:));
%Y(Y<0)=0;




