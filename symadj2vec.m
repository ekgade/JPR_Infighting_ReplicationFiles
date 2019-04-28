function b=symadj2vec(c,idex)
% converts symmetric NxN adjacency matrix to vector of length N(N-1)/2
% of upper diagonal elements (no diagonal)

if nargin<2
    idex=[];
end
a=c;
a(idex,:)=[];
a(:,idex)=[];
N=size(a,1);
n=0;
b=zeros(N*(N-1)/2,1);
for j=1:N
    for k=1:j-1
        n=n+1;
        b(n)=a(k,j);
    end
end
