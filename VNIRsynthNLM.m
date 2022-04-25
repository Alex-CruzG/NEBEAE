function [P,A,Yo]=VNIRsynthNLM(Npixels,ModelType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [Y,P,A,G,u,t]=VNIRsynthNLM(N,Npixels,SNR,PSNR,ModelType)
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
% Npixels --> numbers of pixels in x & y axes
% ModelType --> 0 (LMM-Default), 1 (FM), 2 (GBM), 3 (PNMM) and 4 (MMM)
%
% OUTPUTS
% P --> matrix of end-members 186 x N
% A --> matrix of abundances of N x (Npixels*Npixels)
% Yo --> matrix of measurements
% Yo = P*A+G
%
% April/2023
% DUCD
%
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthetic VNIR Dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=3
Nsamp=Npixels;
x=1:Npixels;
y=1:Npixels;
[xx,yy]=meshgrid(x,y);
K=Npixels*Npixels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End-members extracted from:
% Ines A. Cruz-Guerrero et al, 
% ``Classification of Hyperspectral In Vivo Brain Tissue Based on Linear Unmixing
% Applied Sciences (2020)
% https://doi.org/10.3390/app10165686
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load EndMembersVNIR;
L=size(P,1);
P1=P(:,1); P1=P1/max(P1); 
P2=P(:,2); P2=P2/max(P2); 
P3=P(:,3); P3=P3/max(P3); 


L=256;
n=1:L;
P1=exp(-(n-50).^2 / (2*25^2))+0.3*exp(-(n-100).^2 / (2*20^2))+randn(1,max(n))*0.01;
P2=exp(-(n-108).^2 / (2*70^2))+0.7*exp(-(n-148).^2 / (2*50^2))+randn(1,max(n))*0.01;
P3=0.4*exp(-(n-150).^2 / (2*50^2))+exp(-(n-200).^2 / (2*25^2))+randn(1,max(n))*0.01;
P1=P1'/max(P1);
P2=P2'/max(P2);
P3=P3'/max(P3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Abundance Maps and Nonlinear Interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aa1=rand(Nsamp,Nsamp);
aa2=rand(Nsamp,Nsamp);
aa3=rand(Nsamp,Nsamp);

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

  
Yy=zeros(Nsamp,Nsamp,L);
Gg=zeros(Nsamp,Nsamp,L);

for i=1:Nsamp
    for j=1:Nsamp
        
        if ModelType==1
             g=(a1(i,j)*P1).*(a2(i,j)*P2) + (a1(i,j)*P1).*(a3(i,j)*P3) + (a2(i,j)*P2).*(a3(i,j)*P3);
        elseif ModelType==2
             gamma=rand(3,1);
             g=(a1(i,j)*P1).*(a2(i,j)*P2)*gamma(1)+(a1(i,j)*P1).*(a3(i,j)*P3)*gamma(2)+(a2(i,j)*P2).*(a3(i,j)*P3)*gamma(3);
        elseif ModelType==3
            xi=rand*0.6-0.3;
             g=(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3).*(a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3)*xi;
        else
            g=0;
        end
        y=a1(i,j)*P1 + a2(i,j)*P2 + a3(i,j)*P3; 
                    

        if ModelType==4
            d=0.1+0.2*randn;
            sy=sum(y)
            x=y/sy;
            Yy(i,j,:)=sy*((1-d)*x)./(1-d*x);
        else 
            Yy(i,j,:)=y+g;
        end
   
    end
end

Yo=reshape(Yy,K,L)';
P=[P1 P2 P3];
A=[reshape(a1,1,K);reshape(a2,1,K); reshape(a3,1,K)];








