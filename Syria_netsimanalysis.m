function Syria_netsimanalysis(pars)



inputFile=pars.outFileFull;
load(inputFile,'netall','netav','adjmat','nodelabels','x','idScores')
outputFile=inputFile;



eigOrder=pars.eigOrder;
modmatrix=pars.modmatrix;
identeigvec=pars.identeigvec;
nrun=pars.nrun;
lam=pars.lam;
idNames=pars.idNames;
idNamesFull=pars.idNamesFull;

if size(idNames,2)~=1
    idNames=idNames';
end


K=size(idScores,2);
degs=sum(adjmat,2);
N=length(degs);
if ~isempty(lam)
    M=length(lam);
else
    M=1;
end

corrs=zeros(nrun,K,M);
assorts=zeros(nrun,K,M);
pvalCorr=zeros(K,M);
pvalAssort=zeros(K,M);
rObs=zeros(K,1);
thetaObs=zeros(K,1);
rSimAv=zeros(K,M);
thetaSimAv=zeros(K,M);
ev1=zeros(nrun,K,M);
ev1Av=zeros(K,M);
netMedian=zeros(N,N,M);
pvalStand=zeros(K,1);
rSimDev=zeros(K,M);
thetaSimDev=zeros(K,M);
for kk=1:K
    x=idScores(:,kk);
    [corrObs,pvalStand(kk)]=...
        modularity_corr(adjmat,x,eigOrder,identeigvec,modmatrix);
    rObs(kk)=corrObs;
    [~,assortObs]=assortscalar(adjmat,x);
    thetaObs(kk)=assortObs;
    for mm=1:M
        for jj=1:nrun
            simnet=netall(:,:,jj,mm);
            [corrs(jj,kk,mm),~,~,evdiag]=...
                modularity_corr(simnet,x,eigOrder,identeigvec,modmatrix);
            ev1(jj,kk,mm)=evdiag(1,1);
            [~,assorts(jj,kk,mm)]=assortscalar(simnet,x);
        end
        nsigCorr=2*min(length(find(corrs(:,kk,mm)<corrObs)),...
            length(find(corrs(:,kk,mm)>=corrObs)));
        pvalCorr(kk,mm)=nsigCorr/nrun;
        nsigAssort=2*min(length(find(assorts(:,kk,mm)<assortObs)),...
            length(find(assorts(:,kk,mm)>=assortObs)));
        pvalAssort(kk,mm)=nsigAssort/nrun;
        rSimAv(kk,mm)=mean(corrs(:,kk,mm));
        rSimDev(kk,mm)=std(corrs(:,kk,mm));
        thetaSimAv(kk,mm)=mean(assorts(:,kk,mm));
        thetaSimDev(kk,mm)=std(assorts(:,kk,mm));
        ev1Av(kk,mm)=mean(ev1(:,kk,mm));
        netMedian(:,:,mm)=median(netall(:,:,:,mm),3);
    end
end

if M==1
    summTab=table(idNames,rObs,pvalStand,rSimAv,rSimDev,pvalCorr,...
        thetaObs,thetaSimAv,thetaSimDev,pvalAssort,ev1Av);
else
    summTab=[];
end

netCentral=round(netav);
[~,~,U,S]=modularity_corr(netCentral,x,eigOrder,identeigvec,modmatrix);

outputFile=[outputFile '_sum'];
save(outputFile,'pars','corrs','assorts','pvalCorr','pvalAssort',...
    'rObs','thetaObs','ev1','netav','netMedian','netCentral','pvalStand',...
    'summTab')

if ~isempty(summTab)
    writetable(summTab,[outputFile '.xlsx'])
end

d1=U(:,1);
d2=U(:,2);
d3=U(:,3);


axsign=[1 1];
d1=axsign(1)*d1;
d2=axsign(2)*d2;

if identeigvec==1
    identdim=d1;
elseif identeigvec==2
    identdim=d2;
elseif identeigvec==3
    identdim=d3;
end

ixISIL=find(strcmp('ISIL',nodelabels),1); % place ISIL on positive side
if identdim(ixISIL)<0
    identdim=-identdim;
end

numId=length(idNames);
for idIndex=1:numId
    tp=idScores(:,idIndex);
    idTypeFull=idNamesFull{idIndex};
    idType=idNames{idIndex};
    [cmax,pval]=corr(identdim,tp);


    % calculate assortativity
    [~,assortcoef]=assortscalar(netCentral,tp);
    
    
    
    
end
