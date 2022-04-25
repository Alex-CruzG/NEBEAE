function [E_t,A_t, P_t,X,l]=unmix(Y,p,itmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [E_t,A_t, P_t,X,l]=unmix(Y,p,itmax)
% 
% Blind Nonlinear Unmixing Based in the paper:
% 
% Q. Wei, M. Chen, J. Tourneret and S. Godsill, 
% "Unsupervised Nonlinear Spectral Unmixing Based on a Multilinear Mixing Model," 
% IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 8, pp. 4534-4544, Aug. 2017, 
% doi: 10.1109/TGRS.2017.2693366.
%
% Nicolas Mendoza Chavarria & DUCD
% UASLP
% 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L,N]=size(Y);
%% valores iniciales
Eo=vca(Y,p); %matrix de perfiles
Po=zeros(N,1);%vector de probabilidad 
[Eo,Ao]=LMM(Y,p,Eo);%matriz de abundancias;%matriz de abundancias
%%
E_t=Eo;
P_t=Po;
A_t=Ao;
e=1e-3;
l(1)=fL(E_t,A_t,P_t,Y);
    %disp(['L_o= ', num2str(l(1))]);
for i=2:itmax+1
        A_t1= mina(Y,E_t,A_t,P_t);%update A Algo1
        if (fL(E_t,A_t1,P_t,Y) <= fL(E_t,A_t,P_t,Y) )
            A_t=A_t1;
            %disp(' cambio At');
        end

        P_t= upP(Y, E_t, A_t);%update P
%         figure(1)
%         hist(P_t)

        E_t1=mine(E_t,A_t,Y,P_t); %update E algo2
        if ( fL(E_t1,A_t,P_t,Y) <= fL( E_t,A_t,P_t,Y) )
            E_t=E_t1;
            %disp('cambio Et')
        end
        l(i)=fL(E_t,A_t,P_t,Y);
        %disp(['L= ', num2str(l(i))]);
        if(abs( (l(i)-l(i-1)) / l(i-1) )<e)
%             disp(['error = ',num2str(((l(i)-l(i-1))/l(i-1)))]);
%             disp(['iteracion =',num2str(i-1)]);
            break;
        end
    %disp(['it= ', num2str(i)])
    end
    ylm=E_t*A_t;
    for i=1:N
        X(:,i)=((1-P_t(i)).*ylm(:,i))./(1-(P_t(i).*ylm(:,i)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update P
function P_g =upP(X, E, A)
    %optimaztion respect P
    [L,N]=size(X);
    [s,p]= size(E);
    Y=E*A;%Y=Ea
    %Y=Y./sum(Y);
    P_g=zeros(N,1);
    for i=1:N
        t1=Y(:,i)-( Y(:,i).*X(:,i) );
        den=norm(t1)^2;%||y- y.*x||_2^2
        num=t1'*(Y(:,i)-X(:,i));
    %P_g(i)=projsplx(num/den);%Ec(14)
     P_g(i)=min(1,num/den);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update A
function a_g= mina(Y,Eo,Ao,Po)
    %Algorithm 1 Minimization With Respect to a
    [L,N]=size(Y);
    [s,p]= size(Eo);
    %% Valores iniciales
    Itm=ones(1,p);
    Idm=ones(L,p);
    a_g=zeros(p,N);
    %%
    for i=1:N
        E_g=Eo.*( (1-Po(i)*Idm) + Po(i)*Y(:,i)*Itm);
        L=norm(E_g'*E_g,'Fro');%lipschitz constant
        g_ga=E_g'* (E_g*Ao(:,i)-Y(:,i));
        a_g(:,i)= projsplx(Ao(:,i)-g_ga/L);% gradient projection update
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update E
function E_g = mine(E, A, X, P) 
% Algoritmhm 2
% optimization respect E
% E endmember matrix
% A abundance matrix
% X Hyperspectral image
% P probability
    [p,N]=size(A);
    [L,p]=size(E);

    for i= 1:N
        A_i(:,:,i)=( (1-P(i)) *ones(L,1) +P(i)*X(:,i)) *A(:,i)';% 1
        gEp(:,:,i)=(((E.*A_i(:,:,i))*ones(p,1)-X(:,i))*ones(1,p)).*A_i(:,:,i);% 2 EC(18)
    end

     gE=sum(gEp,3);
     ge=zeros(p,p);
     sn=zeros(p,p,N); 
   for j=1:L
        for i=1:N
            sn(:,:,i)=(A_i(j,:,i).*A_i(j,:,i)');
        end
        ge(:,:)=sum(sn,3);
        fnorm=norm(ge,'fro');
        E_g(j,:)= E(j,:) - ( gE(j,:)/fnorm );
   end
    E_g=min(E_g,ones(L,p));
    E_g=max(E_g,zeros(L,p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% vca
 function [M,indices,snrE]=vca(R,p)
[L,N]=size(R);
SNRth=15+10*log10(p);
rm=mean(R,2);
rzm=R-repmat(rm,1,N);
[U,S,V]=svds(rzm*rzm'./N,p);
rd=U'*rzm;
pr=sum(R(:).^2)/N;
prp=sum(rd(:).^2)/N +rm'*rm;
snrE=abs(10*log10((prp-(p/L)*pr)/(pr-prp)));

if (snrE>SNRth)
    d=p;
    [Ud,S,V]=svds(R*R'./N,d);
    X=Ud'*R;
    u=mean(X,2);
    Y=X ./ repmat( sum( X .* repmat(u,[1 N]) ) ,[d 1]);
else
    d=p-1;
      r_line=mean(R')';
    Ud=pca((R-r_line)*(R-r_line)'./N,'NumComponents',d);
    R_r=R-repmat(r_line,1,N);
    X=Ud'*R_r;
    c=zeros(N,1);
    for j=1:N
        c(j)=norm(X(:,j));
    end
    c=repmat(max(c),1,N);
    Y=[X;c];
end
eu=zeros(p,1);
eu(p)=1;
A=zeros(p,p);
A(:,1)=eu;
I=eye(p);
k=zeros(N,1);
for i=1:p
    w=rand(p,1);
    tmpnum=(I-A*pinv(A))*w;
    f=tmpnum/norm(tmpnum);
    v=f'*Y;
    k=abs(v);
    [unused,k]=max(k);
    A(:,i)=Y(:,k);
    indices(i)=k;
end
if (snrE>SNRth)
    M=Ud*X(:,indices);
else
    M=Ud*X(:,indices)+repmat(r_line,1,p);
end
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L
function l = fL(E,A,P,X)
[L,p]= size(E);
[s,N]= size(A);
y=E*A;
t= X - (1-P').*y - ( (P'.*y) .*X);
d=vecnorm(t).^2;
l=sum(d);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% projeccion
function x = projsplx(y)
% project an n-dim vector y to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
%
% Jan. 14, 2011.
m = length(y); bget = false;
s = sort(y,'descend'); tmpsum = 0;
for ii = 1:m-1
    tmpsum = tmpsum + s(ii);
    tmax = (tmpsum - 1)/ii;
    if tmax >= s(ii+1)
        bget = true;
        break;
    end
end
    
if ~bget, tmax = (tmpsum + s(m) -1)/m; end;
x = max(y-tmax,0);
return;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=FCLS(X,P0)
    [L,N]=size(X);
    p=size(P0,2);
    A=zeros(p,N);
    delta = 1/(10*max(max(P0)));
    for i=1:N
        s = [delta.*X(:,i);1];
        M = [delta.*P0; ones(1,p)];
        A(:,i)= NCLS(s, M, -1e-6);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [abundance]=NCLS(x, MatrixZ, tol)
    % input MatrixZ is the signatures of endmembers. It is of size [bands p].
    % input x is the signature whose abundance is to be estimated.
    % output abundance is the abundance of each material in r1. It is of size [p 1]. % This function is written according to Dr. Chang?s first book , P 47
    M=size(MatrixZ,2);
    R=zeros(M,1);
    P=ones(M,1);
    invMtM=(MatrixZ'*MatrixZ)^(-1);
    Alpha_ls=invMtM*MatrixZ'*x;
    Alpha_ncls=Alpha_ls;
    min_Alpha_ncls=min(Alpha_ncls);
    j=0;
    while(min_Alpha_ncls<-tol && j<500) 
        j = j+1;
        for II=1:M
            if((Alpha_ncls(II)<0)&&(P(II)==1)) R(II)=1;
                P(II)=0;
            end %%% end of if (Alpha_ncls(II)<0) 
        end % end of for II=1:M
        S = R;
    goto_step6=1; 
    counter = 0; 
    while(goto_step6==1)
        index_for_Lamda = find(R==1);
        Alpha_R = Alpha_ls(index_for_Lamda);
        Sai = invMtM(index_for_Lamda,index_for_Lamda);
        inv_Sai = (Sai)^(-1); % remember inversion of Sai 
        Lamda=inv_Sai*Alpha_R;
        [max_Lamda,index_Max_Lamda]=max(Lamda); 
        counter = counter+1;
        if (isempty(max_Lamda))
            break;
        end
        if ( max_Lamda<=0 || counter == 200 )
            break; 
        end
        temp_i = inv_Sai; % simplify the inversion of matrix 
        temp_i(1,:) = inv_Sai(index_Max_Lamda,:);
        if  (index_Max_Lamda>1)
            temp_i(2:index_Max_Lamda,:) = inv_Sai(1:index_Max_Lamda-1,:); 
        end
        inv_Sai_ex = temp_i;
        inv_Sai_ex(:,1) = temp_i(:,index_Max_Lamda);
        if index_Max_Lamda>1
            inv_Sai_ex(:,2:index_Max_Lamda) = temp_i(:,1:index_Max_Lamda-1); 
        end
        inv_Sai_next = inv_Sai_ex(2:end,2:end) - inv_Sai_ex(2:end,1)*inv_Sai_ex (1,2:end)/inv_Sai_ex(1,1);
        P(index_for_Lamda(index_Max_Lamda))=1; 
        R(index_for_Lamda(index_Max_Lamda))=0; 
        index_for_Lamda(index_Max_Lamda) = [];
        Alpha_R = Alpha_ls(index_for_Lamda); 
        Lamda=inv_Sai_next*Alpha_R;
        Phai_column = invMtM(:,index_for_Lamda);
        if (size(Phai_column,2)~=0) 
         Alpha_s=Alpha_ls-Phai_column*Lamda;
        else
            Alpha_s=Alpha_ls;
        end
        goto_step6=0;
        for II=1:M
            if ((S(II)==1)&&(Alpha_s(II)<0))
                P(II)=0;
                R(II)=1;
                goto_step6=1;
            end
        end
    end % end of while (gotostep6==1)
    index_for_Phai = find(R==1);
    Phai_column = invMtM(:,index_for_Phai);
    if (size(Phai_column,2)~=0) 
        Alpha_ncls=Alpha_ls-Phai_column*Lamda;
    else
            Alpha_ncls=Alpha_ls;
    end
    min_Alpha_ncls=min(Alpha_ncls); 
end % end of while
abundance=zeros(M,1); 
for II=1:M
    if (Alpha_ncls(II)>0) 
        abundance(II)=Alpha_ncls(II);
    end
end
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,Ao] = LMM(Y,p,Eo)
[L,N]=size(Y);
% LMM
initcond=5;             % Initial condition of end-members matrix
rho=0.000001;               % Similarity weight in end-members estimation
lambda=0.0;             % Entropy weight for abundance estimation
epsilon=1e-3;
 maxiter=50;
parallel=0;
normalization=1;
downsampling=0.25;       % Downsampling in end-members estimation
disp_iter=0;             % Display partial performance in BEAE          % Display partial performance in BEAE
parameters=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,normalization,disp_iter];
[P,Ao]=BEAE12(Y,p,parameters,Eo,1);
end
