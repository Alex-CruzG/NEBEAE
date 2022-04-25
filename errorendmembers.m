function Ep=errorendmembers(Po,P)

N=size(Po,2);
Error=zeros(N,N);

for i=1:N
    for j=1:N
        Poi=Po(:,i)/sum(Po(:,i));
        Pj=P(:,j)/sum(P(:,j));
        Error(i,j)=norm(Poi-Pj,'Fro')/norm(Poi,'Fro');
    end
end

Ep=min(Error(:));