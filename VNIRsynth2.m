function [Y,P,A,G]=VNIRsynth3(N,Npixels,SNR,PSNR,ModelType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [Y,P,A,G,u,t]=VNIRsynth3(N,Npixels,SNR,PSNR,ModelType)
%
% Synthetic mFLIM dataset with linear or non-linear mixing models:
%
% 1) Linear Mixing Model (LMM)
%
% y_k = \sum_{n=1}^N a_{k,n} p_n + v_k         v_k \sim N(0,\sigma^2 I)
%
% 2) Fan Model (FM) --> Fan et al., 2009
%
% y_k = \sum_{n=1}^N a_{k,n} p_n 
%         + \sum_{n=1}^{N-1} \sum_{m=n+1}^N (a_{k,n}p_n) \odot (a_{k,m}p_m) + v_k
%
% 3) Generalized Bilinear Model (GBM) --> Halimi et al., 2011)
%
% y_k = \sum_{n=1}^N a_{k,n} p_n 
%         + \sum_{n=1}^{N-1} \sum_{m=1}^N \gamma_{n,m} (a_{k,n}p_n) \odot (a_{k,m}p_m) + v_k
%
% 4) Postnonlinear Mixing Model (PNMM) --> Altmann et al., 2012
%
% y_k = \sum_{n=1}^N a_{k,n} p_n 
%         + \sum_{n=1}^{N} \sum_{m=1}^N \xi (a_{k,n}p_n) \odot (a_{k,m}p_m) + v_k
%
% 5) Multilinar Mixing Model (MLM)--> R. Heylen and P. Scheunders., 2016
%  x_k=(1-P_k)*y_k / (1-P_k*y_k) | y_k= sum_{i}a_{i} *w_{i} | w_{i}= p_i/(P_i*p_i +1 - P_i)
%   p_i endmember; P_i probability
%
% INPUTS
% N --> Order of multi-exponential model \in [2,4]
% Npixels --> numbers of pixels in x & y axes
% SNR --> SNR of Gaussian noise (dB)
% PSNR --> PSNR for Shot noise (dB)
% ModelType --> 0 (LMM-Default), 1 (FM), 2 (GBM) and 3 (PNMM)
%
% OUTPUTS
% Y --> matrix of fluorescence decays of size 186 x (Npixels*Npixels)
% P --> matrix of end-members 186 x N
% A --> matrix of abundances of N x (Npixels*Npixels)
% G --> matrix of nonlinear interaction 186 x (Npixels*Npixels) or Prob
% vector  if ModelType==5
% u --> vector of laser input
% t --> vector of time samples
%
% Y = P*A+G
%
% Nov/2020
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
    ModelType=1;
elseif nargin<5
    ModelType=1;
elseif nargin<4
    PSNR=0;
    ModelType=0;
elseif nargin<3
    PSNR=0; SNR=0; ModelType=0;
elseif nargin<2 
    PSNR=0; SNR=0; Npixels=60;ModelType=1;
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
load EndMembersVNIR;
L=size(P,1);
P1=P(:,1); %P1=P1/sum(P1); P1=P1*1000;
P2=P(:,2); %P2=P2/sum(P2); P2=P2*1000;
P3=P(:,3); %P3=P3/sum(P3); P3=P3*1000;
P4=P(:,4); %P4=P4/sum(P4); P4=P4*1000;
maxvalue=max(max([P1 P2 P3 P4]));
%if ModelType==5
    P1=P1/maxvalue;
    P2=P2/maxvalue;
    P3=P3/maxvalue;
    P4=P4/maxvalue;
%end
gamma=[0.95 0.8 0.95 0.5 0.6 0.8]; %Mixing coefficients in GBM
xi=0.3;                              %Scaling coefficient in PNMM
prob=ones(Npixels,Npixels)*0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Abundance Maps and Nonlinear Interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if N==2,

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
    for i=1:Nsamp,
        for l=1:Nsamp,
            a1(i,l)=aa1(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a2(i,l)=aa2(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a3(i,l)=aa3(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
            a4(i,l)=aa4(i,l)/(aa1(i,l)+aa2(i,l)+aa3(i,l)+aa4(i,l));
        end
    end
end
    

Yy=zeros(Nsamp,Nsamp,L);
Gg=zeros(Nsamp,Nsamp,L);

for i=1:Nsamp
    for j=1:Nsamp
        if N==2
            if ModelType==2
                g=(a1(i,j)*P1).*(a2(i,j)*P2);
            elseif ModelType==3,
                g=(a1(i,j)*P1).*(a2(i,j)*P2)*gamma(1);
            elseif ModelType==4,
                g=(a1(i,j)*P1 + a2(i,j)*P2).*(a1(i,j)*P1 + a2(i,j)*P2)*xi;
            else
                g=0;
            end
            y=a1(i,j)*P1 + a2(i,j)*P2;
            
        elseif N==3,
            if ModelType==2,
                 g=(a1(i,j)*P1).*(a2(i,j)*P2) + (a1(i,j)*P1).*(a3(i,j)*P3) + (a2(i,j)*P2).*(a3(i,j)*P3);
            elseif ModelType==3,
                 g=(a1(i,j)*P1).*(a2(i,j)*P2)*gamma(1)+(a1(i,j)*P1).*(a3(i,j)*P3)*gamma(2)+(a2(i,j)*P2).*(a3(i,j)*P3)*gamma(3);
            elseif ModelType==4,
                 g=(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3).*(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3)*xi;
            else
                g=0;
            end
            y=a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3; 
        elseif N==4
            if ModelType==2
                g1=(a1(i,j)*P1).*(a2(i,j)*P2) + (a1(i,j)*P1).*(a3(i,j)*P3) + (a1(i,j)*P1).*(a4(i,j)*P4);
                g=g1 + (a2(i,j)*P2).*(a3(i,j)*P3) + (a2(i,j)*P2).*(a4(i,j)*P4) + (a3(i,j)*P3).*(a4(i,j)*P4);
            elseif ModelType==3
                g1=(a1(i,j)*P1).*(a2(i,j)*P2)*gamma(1) + (a1(i,j)*P1).*(a3(i,j)*P3)*gamma(2) + (a1(i,j)*P1).*(a4(i,j)*P4)*gamma(3);
                g=g1 + (a2(i,j)*P2).*(a3(i,j)*P3)*gamma(4) + (a2(i,j)*P2).*(a4(i,j)*P4)*gamma(5) + (a3(i,j)*P3).*(a4(i,j)*P4)*gamma(6);
            elseif ModelType==4
                g=(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3 + a4(i,j)*P4).*(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3 + a4(i,j)*P4)*xi;
            else
                g=0;
            end
            y=a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3 + a4(i,j)*P4; 
        end
        if ModelType==5
            x=y./sum(y);
            y=((1-prob(i,j))*x)./(1-(prob(i,j)*x));
        end
        Gg(i,j,:)=g;     
        Yy(i,j,:)=y+g;
    end
end

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

%Yy=Yy+Gg;

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
Yo=reshape(Yy,K,L)';
Ym=sum(Yo);
Y=Yo;%./repmat(Ym,[L,1]);
if ModelType==5
Y=Yo;
end
  




