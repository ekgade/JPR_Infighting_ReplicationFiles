function [r,p,U,S]=modularity_corr(adjmat,x,eigOrder,identeig,modmatrix)

if nargin<3
    eigOrder='la';
end
if nargin<4
    identeig=1;
end
if nargin<5
    modmatrix='modularity';
end

degs=sum(adjmat,2);

if isequal(modmatrix,'laplacian')
    degdiag=diag(degs);
    nlap=degdiag^(-1/2)*adjmat*degdiag^(-1/2); % calculate norm. Laplacian
    nlap=(nlap+nlap')/2; % symmetrize to get rid of rounding errors
    [U,S]=eigs(nlap,3,eigOrder);
elseif isequal(modmatrix,'modularity')
    degsum=sum(degs);
    B=adjmat-degs*degs'/degsum;
    [U,S]=eigs(B,3,eigOrder);
end

identdim=U(:,identeig);
[r,p]=corr(identdim,x);

