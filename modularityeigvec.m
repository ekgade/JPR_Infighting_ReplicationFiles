function [eigVec,eigVal] = modularityeigvec(A,numeigs,eigOrder)
%% modularityeigvecs
% Calculates eigenvalues and eigenvectors of symmetric adjacency matrix

if nargin < 2 || isempty(numeigs)
    numeigs=4;
end
if nargin < 3 || isempty(eigOrder)
    eigOrder='la';
end


% check if symmetric
if ~isempty(find(A-A'> 100*eps, 1))
    warndlg('adjacency matrix not symmetric')
end

B=modularitymatrix(A);
B=(B+B')/2; % symmetrize to get rid of roundoff error

[eigVec,S]=eigs(B,numeigs,eigOrder); % sort in specified algebraic order
eigVal=diag(S);

% ensure all components of zero eigenvalue are exactly the same
nzero=find(abs(diag(S))<100*eps);
if ~isempty(nzero)
    eigVec(:,nzero)=repmat(mean(eigVec(:,nzero),1),size(A,1),1);
end



