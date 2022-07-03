function [A_est, time, index]=SVMAX(X,N)     
%=====================================================================
% Programmers: 
% Tsung-Han Chan, E-mail: thchan@ieee.org  
% A. ArulMurugan, E-mail: aareul@ieee.org
% Date: Sept, 2010
%======================================================================
% A implementation of SVMAX
% [A_est time index]=SVMAX(X,N)
%======================================================================
%  Input
%  X is M-by-L data matrix where M is the spectral bands (or observations) and L is the number of pixels (data length).   
%  N is the number of endmembers (or sources).
%----------------------------------------------------------------------
%  Output
%  A_est is M-by-N: estimated endmember signatures (or mixing matrix) obtained by SVMAX.
%  time is the computation time (in secs). 
%  index is the set of indices of the pure pixels identified by SVMAX
%========================================================================

t0 = clock;
[M,L] = size(X);
d = mean(X,2);
U = X-d*ones(1,L);
OPTS.disp = 0;
[C D] = eigs(U*U',N-1,'LM',OPTS);
Xd_t = C'*U;
%=====SVMAX algorithm=========
A_set=[]; Xd = [Xd_t; ones(1,L)]; index = []; P = eye(N);                         
for i=1:N
    [val ind]=max(sum(abs(P*Xd).^2).^(1/2));    
    A_set = [A_set Xd(:,ind)];                            
    P = eye(N) - A_set*pinv(A_set);                       
    index = [index ind];                                        
end
A_est=C*Xd_t(:,index)+d*ones(1,N);
time = etime(clock,t0);

