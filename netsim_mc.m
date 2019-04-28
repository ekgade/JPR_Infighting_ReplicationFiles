function [netav,netall]=netsim_mc(pars,saveOutput)

if nargin<2
    saveOutput=0;
end


outFileFull=pars.outFileFull;

idNames=pars.idNames;
idNamesFull=pars.idNamesFull;
idScores=pars.idScores;
tgrps=pars.tgrps;

lam=pars.lam;
nrun=pars.nrun;
idIdx=pars.idIdx;
eigOrder=pars.eigOrder;
binarize=pars.binarize;
modmatrix=pars.modmatrix;
binNetForPlot=pars.binNetForPlot;
identeigvec=pars.identeigvec;

adjmat=pars.net;
nodelabels=pars.netgroups;
nettype=pars.nettype;

degs=sum(adjmat,2);
izdeg=find(degs==0);
if ~isempty(izdeg)
    degs(izdeg)=[];
    nodelabels(izdeg)=[];
    adjmat(izdeg,:)=[];
    adjmat(:,izdeg)=[];
    idScores(izdeg,:)=[];
end

x=idScores(:,idIdx);

% get observed correlations, assortativity
[corrObs,pObs]=modularity_corr(adjmat,x,eigOrder,identeigvec,modmatrix);
[~,assortObs]=assortscalar(adjmat,x);

        
N=length(degs);
if ~isempty(lam)
    M=length(lam);
else
    M=1;
end
netinp=zeros(N);
corrs=zeros(nrun,M);
assorts=zeros(nrun,M);
netall=zeros(N,N,nrun,M);
netav=zeros(N,N,M);
pvalCorr=zeros(M,1);
pvalAssort=zeros(M,1);
netExpect=zeros(N,N,M);
for mm=1:M
    nettot=zeros(N);
    if isempty(lam)
        lamval=[];
    else
        lamval=lam(mm);
    end
    jj=1;
    while jj<=nrun
        if pars.simtype==0 || pars.simtype==1
            [netevol,netExp]=jonetsim(degs,x,lamval,netinp);  % homphily
        elseif pars.simtype==-1
            [netevol,netExp]=jonetsim_anti(degs,x,lamval,netinp); % heterophily
        end
        netall(:,:,jj,mm)=netevol;
        [r,p,U,S]=modularity_corr(netevol,x,eigOrder,identeigvec,modmatrix);
        corrs(jj,mm)=r;
        [~,assort]=assortscalar(netevol,x);
        assorts(jj,mm)=assort;
        nettot=nettot+netevol;
        jj=jj+1;
    end
    netav(:,:,mm)=nettot/nrun;
    nsigCorr=length(find(abs(corrs(:,mm))>=abs(corrObs)));
    pvalCorr(mm)=nsigCorr/nrun;
    nsigAssort=2*min(length(find(assorts(:,mm)<assortObs)),...
        length(find(assorts(:,mm)>=assortObs)));
    pvalAssort(mm)=nsigAssort/nrun;
    netExpect(:,:,mm)=netExp;
end

if saveOutput
    save(outFileFull,'pars','netall','netav','pvalCorr','pvalAssort',...
        'corrs','assorts','corrObs','netExpect','adjmat','assortObs',...
        'nodelabels','idScores','x')
end
