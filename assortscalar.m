function [r,rcoef] = assortscalar(A,X)
%% assortscalar
% Calculates network assortativity for scalar node attribute

%% Code

A=(A+A')/2; % symmetrize A

degs=sum(A,2);
degsum=sum(degs);
B=A-degs*degs'/degsum;
r=X'*B*X/degsum;

BPerfect=diag(degs)-degs*degs'/degsum;
rPerfect=X'*BPerfect*X/degsum;
rcoef=r./rPerfect;

