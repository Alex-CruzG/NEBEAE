%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Unmixing with Nonlinear and Linear Models for in-vivo barin VNIR Dataset --> Figure 8
%
% ``Nonlinear Extended Blind End-member and Abundance Extraction for
% Hyperspectral Images''
%  Campos-Delgado D.U. et al, Submitted to Signal Processing (Elsevier)
%
%
% DUCD
% April/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for dataset=1:2

if dataset==1      
    file1='VNIRimagesOp15C1';
    titleL='(a)Op15C1';
    file2='P3GoldenReference.mat';
    file3='Op15C1_Y.jpg';
elseif dataset==2
    file1='VNIRimagesOp20C1';
    titleL='(a)Op20C1';
    file2='P4GoldenReference.mat';
    file3='Op20C1_Y.jpg';
end

load(file1);
load(file2);
 if dataset == 6
     goldenStandardMap=gtMap;
 end
disp('%%%%%%%%%%%%%%%%%%%');
disp(file1);
disp('Dataset loaded!!');
disp('%%%%%%%%%%%%%%%%%%%');
rawImager=rawImage(1:(end-1),1:(end-1),1:end);


%%%%%%%%%%%%%%%%%%%%%
% EBEAE Unmixing
%%%%%%%%%%%%%%%%%%%%%
%%
epsilon=1e-3;            % Threshold for convergence
maxiter=50;              % Maximum number of iteration in alternated least squares approach
parallel=0;              % Parallelization in abundance estimation process
downsampling=0.0;
Image=preProcessedImage; 
fileImage='preProcessedImage';
[Ny,Nx,Nz]=size(Image);
M=reshape(goldenStandardMap,Nx*Ny,1);
Z=reshape(Image,Nx*Ny,Nz)';
I1=find(M==1);
I2=find(M==2);
I3=find(M==3);
I4=find(M==4);
Z1=Z(:,I1);
Z2=Z(:,I2);
Z3=Z(:,I3);
Z4=Z(:,I4);
P1l=mean(Z1,2); 
P2l=mean(Z2,2); 

%%

[L,K]=size(Z);
Y=Z./repmat(sum(Z,1),L,1);
paramvec=[6,0.75,0,epsilon,maxiter,downsampling,0,1,0];
[P3l,~,~,~]=BEAE12(Z3,2,paramvec);
paramvec=[6,0.01,0,epsilon,maxiter,downsampling,0,1,0];
[P4l,~,~,~]=BEAE12(Z4,4,paramvec);
Ppl=[P1l P2l P3l P4l];
Pol=Ppl./repmat(sum(Ppl,1),L,1);

nl=size(Pol,2);
rho=0.1;                   
lambda=0.2;  
initcond=6;
epsilon=1e-2;            % Threshold for convergence
maxiter=50;              % Maximum number of iteration in alternated least squares approach
parallel=0;              % Parallelization in abundance estimation process
downsampling=0.0;
normalization=1;         % Normalization of end-members (sum-to-one)
disp_iter=1;          % Display results of iterative optimization
oae=1;

disp('%%%%%%%%%%%%%%%%%%%%');
disp(fileImage);

paramvec2=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,normalization,disp_iter];
[Pl,al,Sl,Yhl]=BEAE12(Y,nl,paramvec2,Pol,oae);

%%%%%%%%%%%%%%%%%%%%%
% NEBEAE Unmixing
%%%%%%%%%%%%%%%%%%%%%

%%

tic;
paramvec=[6,0.75,0,epsilon,maxiter,downsampling,0,0];
[P1,~,D1,~]=NEBEAE(Z1,2,paramvec);
paramvec=[6,0.75,0,epsilon,maxiter,downsampling,0,0];
[P2,~,D2,~]=NEBEAE(Z2,2,paramvec);
paramvec=[6,0.75,0,epsilon,maxiter,downsampling,0,0];
[P3,~,D3,~]=NEBEAE(Z3,2,paramvec);
paramvec=[6,0.01,0,epsilon,maxiter,downsampling,0,0];
[P4,~,D4,~]=NEBEAE(Z4,4,paramvec);
Pp=[P1 P2 P3 P4];
Po=Pp./repmat(sum(Pp,1),L,1);
n=size(Po,2);

paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,disp_iter];
[P,a,D,S,Yh]=NEBEAE(Y,n,paramvec,Po,oae);

n11=1;
n12=n11+size(P1,2)-1;
n21=n12+1;
n22=n21+size(P2,2)-1;
n31=n22+1;
n32=n31+size(P3,2)-1;
n41=n32+1;
n42=n41+size(P4,2)-1;
m11=1;
m12=m11+size(P1l,2)-1;
m21=m12+1;
m22=m21+size(P2l,2)-1;
m31=m22+1;
m32=m31+size(P3l,2)-1;
m41=m32+1;
m42=m41+size(P4l,2)-1;
a=[sum(a(n11:n12,:),1); sum(a(n21:n22,:),1); sum(a(n31:n32,:),1); sum(a(n41:n42,:),1) ];
al=[sum(al(m11:m12,:),1); sum(al(m21:m22,:),1); sum(al(m31:m32,:),1); sum(al(m41:m42,:),1) ];

%%
n=4;
Aan=zeros(n,K);
for k=1:K
      abunclass=zeros(n,1);
      [Amax,Imax]=max(al(:,k)); abunclass(Imax)=1;
      Aan(:,k)=abunclass;
end

Ac=zeros(Ny,Nx); 
index=1;
for i=1:n,       
        Aa=zeros(Nx*Ny,1);
        eval(['Aa=i*Aan(' num2str(i) ',:);']);
        eval(['AAn' num2str(i) '=reshape(Aa,Ny,Nx);']);
        eval(['Ac=Ac+AAn' num2str(i) ';']);
end;


se3 = strel('disk',1);
Acl=Ac;
Ac1l=imclose(Ac,se3);
Accl=imopen(Ac1l,se3);


Aan=zeros(n,K);
for k=1:K
      abunclass=zeros(n,1);
      [Amax,Imax]=max(a(:,k)); abunclass(Imax)=1;
      Aan(:,k)=abunclass;
end

Ac=zeros(Ny,Nx); 
index=1;
for i=1:n      
        Aa=zeros(Nx*Ny,1);
        eval(['Aa=i*Aan(' num2str(i) ',:);']);
        eval(['AAn' num2str(i) '=reshape(Aa,Ny,Nx);']);
        eval(['Ac=Ac+AAn' num2str(i) ';']);
end

se3 = strel('disk',1);
Ac1=imclose(Ac,se3);
Acc=imopen(Ac1,se3);


subplot(2,5,(dataset-1)*5+1);
imagesc(imread(file3));
title('Approximated RGB Image');
subplot(2,5,(dataset-1)*5+2);
imagesc(generateColorMap(goldenStandardMap));
title('Ground-Truth Map');
subplot(2,5,(dataset-1)*5+3);
imagesc(generateColorMap(Accl)); 
title({titleL,'Classified by EBEAE'});
subplot(2,5,(dataset-1)*5+4);
imagesc(generateColorMap(Acc)); 
title('Classified by NEBEAE');
subplot(2,5,(dataset-1)*5+5);
imagesc(reshape(D,Ny,Nx))
colormap jet
title({'Nonlinear Interaction','Level'});

end

subplot(2,5,5);
colorbar;

