function Ea=errorabundances(Ao,A)

N=size(Ao,1);
Error=zeros(N,N);

for i=1:N
    for j=1:N
        Aoi=Ao(i,:);%/sum(Ao(i,:));
        Aj=A(j,:);%/sum(A(j,:));
        Error(i,j)=norm(Aoi-Aj,'Fro')/norm(Aoi,'Fro');
    end
end

Ea=min(Error(:));