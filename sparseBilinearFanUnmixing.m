function [A,B,S,G,hS]=sparseBilinearFanUnmixing(Y,A,S,maxG,delta,q,hS)
% Model: Y=SA^T+ZB^T
% Method solves: argmin_{S,A,Z,B} 1/2||Y-SA^T-ZB^T||^2_F + h_S ||S||_{q,1}
% s.t. A>0, S>0, 0<G<maxG,
% Y data (PxM) P=#pixels, M=#spectral bands
% A endmembers
% S linear abundances
% G Gamma - 0 < g < maxG
% B bilinear endmembers =[(A_1 * A_2) (A_1 * A_3) ... ]
% Z bilinear abundances Z=[ (G_1 * S_1 * S_2) (G_2 * S_1 * S_3) ... ]
    %Px=50;

%     tic
    [P,~]=size(Y);
    [M,r]=size(A);
    N=r*(r-1)/2;
    maxIters=10000;
    stopCriteria=10^3;%10^-5;
    mainStopCriteria=10^-4;
    alpha=0.5;

    if maxG==0
        G=zeros(P,N);
    else
        G=maxG*rand(P,N);
        G(createB(S)==0)=0;
    end
    
    B=createB(A);
    Z=createB(S);
    Nq=createNq(r);
    Aold=A;Sold=S;Gold=G;
    %Z=createZ(G,S);
    %hS=repmat(hS,1,r);
    
    u=0.01*norm(Y,'fro').^2/(numel(Y));
    
	Ad=[A;delta*ones(1,r)];
    Bd=[B;zeros*ones(1,N)];
    Yd=[Y delta*ones(P,1)];
    global checkConv
    costIter=1;
    cols=[ 'm' 'c' 'r' 'b' 'k']; str=num2str(cputime); s =str2num(str(end))+1;c=ceil(s/2); plotcol=cols(c);
    if isempty(hS)
        [hS, all_lambdas, all_ebic]=estimateLambdaWithEbicFan(Yd,Ad,Bd,S,maxG,q,alpha,0.02);
    end
    for iter=1:maxIters
        checkConv=5;
        if maxG>0
            G=sparseBilinearUnmixingFanGStepSparsePen(Y,A,G,S,Z,stopCriteria,checkConv,maxG,u,q,hS,delta);
        end
        [S,Z]=sparseBilinearUnmixingFanSStep(Yd,Ad,Bd,G,S,hS,q,delta,stopCriteria,checkConv);
        [A,B]=Astep_dyadic(Y,B,G,S,A,stopCriteria/10,Nq,0,checkConv);
        Ad=[A;delta*ones(1,r)];
        Bd=[B;zeros*ones(1,N)];
        if rem(iter,10)==0
            normA=norm(A-Aold,'fro')/norm(A,'fro');
            normS=norm(S-Sold,'fro')/norm(S,'fro');
            normG=0.01*norm(G-Gold,'fro')/max(10^-6,norm(G,'fro'));
            %[normE normA normB]
            if(normA<mainStopCriteria && normS<mainStopCriteria && normG<mainStopCriteria )
                break;
            end
        end
        Aold=A;Sold=S;Gold=G;
        
    end
%     toc
    
end





function [A,B]=Astep_mm(Y,B,S,G,Z,A,stopCriteria,checkConv)
    GZ=G.*Z;
    for i=1:10000
        Aold=A;
        X=max(0,Y-GZ*B');
        A=(A.*(X'*S))./(A*(S'*S));
        B=createB(A);
        if rem(i,checkConv)==0
            normA=norm(A-Aold,'fro')/norm(A,'fro');
            if(normA<stopCriteria)
                return;
            end
        end
    end
end


 function [A,B]=Astep_dyadic(Y,B,G,S,A,stopCriteria,Nq,hA,checkConv)
    Z=G.*createB(S);    
    [M,r]=size(A);
    indx=zeros(r,r-1);
    for i=1:r
        tmp=1:r;
        tmp(i)=[];
        indx(i,:)=tmp;
    end
    N=r*(r-1)/2;
    aN=1:N;
    for i=1:10000
        Aold=A;
        for q=randperm(r)
            nq=Nq(q,:);
            neq=aN;neq(nq)=[];
            mj=indx(q,:);
            Xq=[S(:,mj) Z(:,neq)]*[A(:,mj) B(:,neq)]';
            Yq=Y-Xq;
           
            Zq=[S(:,q) Z(:,nq)]*[ones(M,1) A(:,mj)]';
            aq=calcAq(Yq,Zq);
            %aq=diag(Yq'*Zq)./sum(Zq.^2)';
            aq=max(0,aq);
            A(:,q)=aq; 
            B=updateB(A,B,q);
        end
         if rem(i,checkConv)==0
            normA=norm(A-Aold,'fro')/norm(A,'fro');
            if(normA<stopCriteria)
                return;
            end
            
        end
    end
 end



function aq=calcAq(YmX,Z)
    [~,M]=size(YmX);
    aq=zeros(M,1);
    for m=1:M
        z=Z(:,m);
        aq(m)=YmX(:,m)'*z/(z'*z);
    end
end
 
 
 
function aq=calcSlowEqhA(YmX,Z,hA)
    [~,M]=size(YmX);
    aq=zeros(M,1);
    for m=1:M
        z=Z(:,m);
        aq(m)=YmX(:,m)'*z/(z'*z+hA);
    end
end






% Row q hold indexes in M that have either eq*ei or ei*eq
function Nq=createNq(r)
    N=r*(r-1)/2;
    Nq = tril(ones(r-1));
    Nq(Nq==1) = 1:N;
    Lq=zeros(r);
    Lq(2:end,1:end-1)=Nq;
    Uq=Lq';
    Nq = Lq(:, 1:end-1) + Uq(:, 2:end);
end

function [lambda, all_lambdas, all_ebic]=estimateLambdaWithEbicFan(Y,Ad,Bd,S,maxG,q,alpha,d)
%(Y,Ad,Bd,S,hS,hZ,q,stopCriteria,checkConv)
% Estimate lambda with EBIC
done=0;
%lambdas=[0.001 0.005 0.01 0.015];
%lambdas=[0.0010    0.0023    0.0037    0.0050]*3;

Yo=Y(:,1:end-1);
Ao=Ad(1:end-1,:);


%d=0.02;
lambdas=[0.0010 d 2*d 3*d];

ebic_iter=10000;
checkConv=10;
stopCriteria=10^-4;
[M,r]=size(Ad);
%A0=[Ad Bd];
delta=10;
Z=createB(S);
%u=10*norm(Y)/(numel(Y));
u=0.01*norm(Y,'fro').^2/(numel(Y));
G=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,maxG*rand(size(Z)),S,Z,stopCriteria,checkConv,maxG,u,q,0,delta);

for i=1:4
    S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
    Z=createB(S0{i});
    G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
    ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
end
[a,idx]=min(ebic);
iter=0;
all_lambdas=[lambdas];all_ebic=[ebic];
while done==0
    iter=iter+1;
    if idx==1
        lambdas(1)=max(lambdas(1)-d,0);
        lambdas(4)=lambdas(2); ebic(4)=ebic(2);
        d=(lambdas(4)-lambdas(1))/3;
        lambdas(2)=lambdas(1)+d;
        lambdas(3)=lambdas(1)+2*d;
        for i=1:3
            S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
            Z=createB(S0{i});
            G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
            ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
        end
    end
    if idx==2
        lambdas(4)=lambdas(3); ebic(4)=ebic(3);
        d=(lambdas(4)-lambdas(1))/3;
        lambdas(2)=lambdas(1)+d;
        lambdas(3)=lambdas(1)+2*d;
        for i=[2 3]
            S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
            Z=createB(S0{i});
            G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
            ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
        end
    end
    if idx==3
        lambdas(1)=lambdas(2);ebic(1)=ebic(2);S0{1}=S0{2};
        d=(lambdas(4)-lambdas(1))/3;
        lambdas(2)=lambdas(1)+d;
        lambdas(3)=lambdas(1)+2*d;
        for i=[2 3]
            S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
            Z=createB(S0{i});
            G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
            ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
        end
    end
    if idx==4
        lambdas(1)=lambdas(3);ebic(1)=ebic(3);S0{1}=S0{3};
        lambdas(2)=lambdas(4);ebic(2)=ebic(4);
        lambdas(3)=lambdas(1)+2*d;
        lambdas(4)=lambdas(1)+3*d;
        for i=3:4
            S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
            Z=createB(S0{i});
            G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
            ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
        end
    
        
    end
    [a,idx]=min(ebic);
    
    all_lambdas=[all_lambdas lambdas];
    all_ebic=[all_ebic ebic];
    [a,b]=sort(all_lambdas);
%    figure(1);hold off;plot(all_lambdas(b),all_ebic(b),'x-'); hold on;plot(lambdas,ebic,'rx-');hold off;shg;drawnow
    if abs((lambdas(1)-lambdas(4))/lambdas(4)) < 0.1 || iter>5
        done=1;
        lambda=lambdas(idx);
    end
    %A=A0{1};S=S0{1};
end
[a,b]=min(all_ebic);
lambda=all_lambdas(b)

end

function G=sparseBilinearUnmixingFanGStepSparsePen(Y,A,G,S,Z,stopCriteria,checkConv,maxG,u,q,hS,delta)
    %stopCriteria=10^-5;
    [M,r]=size(A);
    G(Z==0)=0;
    [P,M]=size(Y);
    %nonZeros=logical(ones(size(Z)));
    nonZeros=(Z>0);
    B=createB(A);
    N=r*(r-1)/2;
    indx=zeros(N,N-1);
    %u=0.01;
  	for n=1:N
        tmp=1:N;
        tmp(n)=[];
        indx(n,:)=tmp;
    end
    %SAT=S*A';
    for n=1:N
        nonZ=nonZeros(:,n);
        Pk{n}=Z(nonZ,n)*B(:,n)';
        %[a,b]=hist(Pk{n}(:),100);
        BTB=B(:,n)'*B(:,n);
        tmp=Z(nonZ,n).^2*BTB+u;
        PkPk{n}=1./tmp;
    end
    
    YmSAT=Y-S*A';
    costIter=1;
    for Giter=1:10000
        Gold=G;
        for n=1:N
            nonZ=nonZeros(:,n);
            mn=indx(n,:);
            % Create Rk
            Rk=(G(nonZ,mn).*Z(nonZ,mn))*B(:,mn)';
            % Create Xk
            %Xk=Y(nonZ,:)-SAT(nonZ,:)-Rk;
            Xk=YmSAT(nonZ,:)-Rk;
            g=sum(Pk{n}.*Xk,2).*PkPk{n};
            g=max(0,g);
            G(nonZ,n)=min(maxG,g);
            %cost(costIter)=calcBilFanCost(Y,S,A,G,q,hS,delta,u);costIter=costIter+1;
        end
        %if max(diff(cost))>0 || any(isnan(cost))
        %    foo=0;
        %end
        if rem(Giter,checkConv)==0
            ndG=norm(G-Gold,'fro')/norm(Gold,'fro');
            if(ndG<stopCriteria)
                return;
            end
        end
    end
end

function [S,Z]=sparseBilinearUnmixingFanSStep(Y,A,B,G,S,hS,q,delta,stopCriteria,checkConv)

    [P,r]=size(S);
    %B=createB(A);
    N=r*(r-1)/2;
    Nk=createNq(r);
    aN=1:N;
  	for i=1:r
        tmp=1:r;
        tmp(i)=[];
        indx(i,:)=tmp;
    end
    nonZeros=S>0;
    pS=sum(nonZeros);
    
    hSq=(hS*q);
    qmI=q-1;
    %costIter=1;
    %cost(costIter)=calcBilFanCost(Y(:,1:end-1),S,A(1:end-1,:),G,q,hS,delta);costIter=costIter+1;
    Z=createB(S);
    %SAT=S*A';
    %multSum=@(a,b) sum(a.*b,2);
    for Siter=1:10000
        Sold=S;
        for k=1:r
            nonZ=nonZeros(:,k);
            nk=Nk(k,:);
            nek=aN;nek(nk)=[];
            mk=indx(k,:);
            
            Sn=S(nonZ,:);
            Gn=G(nonZ,:);

            % Create Rk
            %Rk=[S(:,mk) G(:,nek).*Z(:,nek)]*[A(:,mk) B(:,nek)]';
            Rk=[Sn(:,mk) Gn(:,nek).*Z(nonZ,nek)]*[A(:,mk) B(:,nek)]';
            % Create Pk
            %Pk=[G(:,nk).*S(:,mk) ones(P,1)]*[B(:,nk) A(:,k)]';
            Pk=[Gn(:,nk).*Sn(:,mk) ones(pS(k),1)]*[B(:,nk) A(:,k)]';

            %cost(costIter)=calcBilFanCost(Y(:,1:end-1),S,A(1:end-1,:),G,q,hS,delta);costIter=costIter+1;
            s=S(nonZ,k);
            if hS>0
                s=s.*sum(Pk.*(Y(nonZ,:)-Rk),2)./( s.*sum(Pk.*Pk,2)+hSq*exp(qmI*log(s)));
                %s=s.*sum(Pk.*(Y(nonZ,:)-Rk),2)./( s.*sum(Pk.*Pk,2)+hSq*s.^qmI);
            else
                %s=sum(Pk.*(Y-Rk),2)./sum(Pk.*Pk,2);
                s=sum(Pk.*(Y(nonZ,:)-Rk),2)./sum(Pk.*Pk,2);
            end
            %S(:,k)=max(0,s);
            S(nonZ,k)=max(0,s);
            Z=updateB(S,Z,k);
            %cost(costIter)=calcBilFanCost(Y(:,1:end-1),S,A(1:end-1,:),G,q,hS,delta);costIter=costIter+1;

            %for p=1:P
            %    s=S(p,k);
            %    s2(p)=(s*Pk(p,:)*X(p,:)')/(s*Pk(p,:)*Pk(p,:)'+(hS/q)*s^(q-1));
            %end
        end
        if rem(Siter,checkConv)==0
            ndS=norm(S-Sold,'fro')/norm(Sold,'fro');
            if(ndS<stopCriteria/100)
                return;
            end
        end
    end
end

function [B,I,J]=createB(A)
% Creates bilinear endmember combos using Mxr endmember matrix 
% A, where M are the number spectral bands and r is the number of endm.
    [M,r]=size(A);
    N=r*(r-1)/2;
    B=zeros(M,N);
    I=zeros(1,N);
    J=I;
    idx=0;
    for i=1:r-1
        for j=i+1:r
            idx=idx+1;
            B(:,idx)=A(:,i).*A(:,j);
            I(idx)=i;J(idx)=j;
        end
    end
end

function EBIC=EBICfan(Y,S,A,G,alpha)
%
%
B=createB(A);
Z=createB(S);

[T,~]=size(Y);
[M,r]=size(A);
sigma2=1/(T*M)*norm(Y-S*A'-(G.*Z)*B','fro')^2;
de=sum(S(:)>0)+M*r-r^2;
EBIC=M*log(sigma2)+M+1/T*(log(T)+4*alpha*log(M))*de;
end

function [B,I,J]=updateB(A,B,k)
% Creates bilinear endmember combos using Mxr endmember matrix 
% A, where M are the number spectral bands and r is the number of endm.
    [~,r]=size(A);
    N=r*(r-1)/2;
    I=zeros(1,N);
    J=I;
    idx=0;
    for i=1:r-1
        for j=i+1:r
            idx=idx+1;
            if j==k || i==k
                B(:,idx)=A(:,i).*A(:,j);
                I(idx)=i;J(idx)=j;
            end
        end
    end
end
