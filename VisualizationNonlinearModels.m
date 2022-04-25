%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Synthetic 3D Manifold of Nonlinear Mixing Models --> Figure 1
%
% DUCD
% April/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsamples=60;

for index=0:4
    
ModelType=index;

[Po,Ao,Yo]=VNIRsynthNLM(Nsamples,ModelType);
[L,K]=size(Z);
Y=Z;
k=1:L;

if index<3
    subplot(2,3,index+1);
else
    subplot(2,2,index);
end
[U,S,V]=svd(Yo,'econ');
Up=U(:,1:3);
Zp=Up'*Yo;
Pp=Up'*Po;
plot3(Zp(1,:),Zp(2,:),Zp(3,:),'.')
hold on;
for i=1:n
    plot3(Pp(1,i),Pp(2,i),Pp(3,i),'dr');
    text(double(1.005*Pp(1,i)),double(1.005*Pp(2,i)),double(1.005*Pp(3,i)),['EM_' num2str(i)],'FontSize',12)
end
for i=1:(n-1)
    for j=(i+1):n
        plot3([Pp(1,i) Pp(1,j) ],[Pp(2,i) Pp(2,j)],[Pp(3,i) Pp(3,j)],'r-', 'LineWidth',2);
    end
end
view(2); view(3);
view(-20,5)
grid on
if index==0
    title('LMM');
elseif index==1
    title('Fan et al. Model');
elseif index==2
    title('GBM');
elseif index==3
    title('PPNM');
else
    title('MMM');
end
end