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
