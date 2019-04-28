function [a,aExp]=jonetsim_anti(D,x,lam,netinp)

% JONETSIM_ANTI 
% USAGE: a=jonetsim_anti(D,x,lam,netinp)
% INPUTS:
% D: specified node degrees
% x: node scalar attribute
% lam: standard deviation of interaction function
% netinp: initial network (symmetric)
% OUTPUT: 
% a: joint operations network

N=length(D);
M=0.5*sum(D);
if nargin<4
    netinp=zeros(N);
end



if size(D,2)~=1
    D=D';
end
degnow=zeros(N,M+1);
a=zeros(N);

% calculate attribute interaction function values
if ~isempty(lam)
    intFun=zeros(N);
    for j=1:N
        intFun(j,:)=2-exp(-0.5*(x-x(j)).^2/lam^2);
    end
    xInt=intFun-diag(diag(intFun));
else
    xInt=ones(N)-eye(N);
end



for k=1:M
    poss=find(D-degnow(:,k)>0); % nodes with operations remaining
    Nrem=length(poss);
    if Nrem~=1
        opsrem=D-degnow(:,k);
        degProd=opsrem*opsrem';  % degree-degree product
        % interaction probability due to attribute
        probX=degProd.*xInt;
        % normalize and remove upper triangle
        pjoint=tril(probX,-1)/sum(sum(tril(probX,-1)));
        % transform pair probabilities to segments on unit interval
        sep=0;
        ipair=zeros(Nrem*(Nrem-1)/2,2);
        pseg=zeros(Nrem*(Nrem-1)/2,1);
        nloop=0;
        for ii=1:Nrem
            for jj=1:ii-1
                nloop=nloop+1;
                ipair(nloop,:)=[poss(ii) poss(jj)];            
                pseg(nloop)=sep+pjoint(poss(ii),poss(jj));
                sep=pseg(nloop);
            end
        end
        mm=find(pseg-rand>=0,1); % find chosen pair index
        iout=ipair(mm,1);
        jin=ipair(mm,2);
        a(iout,jin)=a(iout,jin)+1; % increment current network
        a(jin,iout)=a(iout,jin);
        degnow(iout,k+1:M+1)=degnow(iout,k)+1; % increment current degree
        degnow(jin,k+1:M+1)=degnow(jin,k)+1;
    end
end

% calculate expected network
degProd=D*D';  
probX=degProd.*xInt;
pjoint=probX/sum(sum(probX));
aExp=2*M*pjoint;

    