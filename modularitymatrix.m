function B = modularitymatrix(A)
%% modularitymatrix
%   Calculates modularity matrix from symmetric adjacency matrix

degs=sum(A,2);
degsum=sum(degs);
B=A-degs*degs'/degsum;



end

