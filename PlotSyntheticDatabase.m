%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Plot of Synthetic Databases --> Figures 2 & 3
%
% ``Nonlinear Extended Blind End-member and Abundance Extraction for
% Hyperspectral Images''
%  Campos-Delgado D.U. et al, Submitted to Signal Processing (Elsevier)
%
%
% DUCD
% April/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SNR=20;     
PSNR=20;
Nsamples=60;
n=4;
ModelType=5;

[Z,Po,Ao,Go]=VNIRsynth(n,Nsamples,SNR,PSNR,ModelType);
[L,K]=size(Z);
k=1:L;

figure;
subplot(2,2,1);
imagesc(reshape(Ao(1,:),Nsamples,Nsamples));
title('(a) Abundance End-member 1');
subplot(2,2,2);
imagesc(reshape(Ao(2,:),Nsamples,Nsamples));
title('(b) Abundance End-member 2');
subplot(2,2,3);
imagesc(reshape(Ao(3,:),Nsamples,Nsamples));
title('(c) Abundance End-member 3');
subplot(2,2,4);
imagesc(reshape(Ao(4,:),Nsamples,Nsamples));
title('(c) Abundance End-member 4');

BloediteWLf;BloediteEMf;
WL2=BloediteWL(BloediteEM>0)
figure;
subplot(2,1,1);
plot(WL2(1:size(Po,1)),Po,'LineWidth',2); grid on;
legend('End-member 1','End-member 2','End-member 3', 'End-member 4');
title('(a) Spectral Library Version 7');
xlabel('Wavelength (\mum)');
axis([min(WL2(1:size(Po,1))) max(WL2(1:size(Po,1))) 0 1]);
size(Po,1)
min(WL2(1:size(Po,1))), max(WL2(1:size(Po,1)))

[Z,Po,Ao,Go]=VNIRsynth2(n,Nsamples,SNR,PSNR,ModelType);
[L,K]=size(Z);
k=1:L;
subplot(2,1,2);
w=linspace(450,950,size(Po,1));
plot(w,Po,'LineWidth',2); grid on;
legend('End-member 1','End-member 2','End-member 3', 'End-member 4');
title('(b) Spectral Responses of Craneotomy Tissue');
xlabel('Wavelength (nm)');
