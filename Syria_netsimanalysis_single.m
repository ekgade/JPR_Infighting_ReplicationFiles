function Syria_netsimanalysis_single(pars)

inputFile=pars.outFileFull;

load(inputFile,'netall','netav','adjmat','nodelabels','x','corrObs',...
    'assortObs')
outputFile=inputFile;

eigOrder=pars.eigOrder;
modmatrix=pars.modmatrix;
binNetForPlot=pars.binNetForPlot;
identeigvec=pars.identeigvec;
nrun=pars.nrun;
lam=pars.lam;
idNames=pars.idNames;
idNamesFull=pars.idNamesFull;
scalName=idNames{pars.idIdx};
idIdx=pars.idIdx;

if size(idNames,2)~=1
    idNames=idNames';
end

ixISIL=find(strcmp('ISIL',nodelabels),1); % used to flip correlation sign


degs=sum(adjmat,2);
N=length(degs);
if ~isempty(lam)
    M=length(lam);
else
    M=1;
end

corrs=zeros(nrun,M);
assorts=zeros(nrun,M);
pvalCorr=zeros(M,1);
pvalAssort=zeros(M,1);
rSimAv=zeros(M,1);
thetaSimAv=zeros(M,1);
ev1=zeros(nrun,M);
ev1Av=zeros(M,1);
netMedian=zeros(N,N,M);
rSimDev=zeros(M,1);
thetaSimDev=zeros(M,1);
lsqerr=zeros(nrun,M);
[rObs,pvalStand]=...
    modularity_corr(adjmat,x,eigOrder,identeigvec,modmatrix);
[~,thetaObs]=assortscalar(adjmat,x);
obsvec=symadj2vec(adjmat);
for mm=1:M
    for jj=1:nrun
        simnet=netall(:,:,jj,mm);
        [corrs(jj,mm),~,Ut,evdiag]=...
            modularity_corr(simnet,x,eigOrder,identeigvec,modmatrix);
        if Ut(ixISIL,1)<0  % flip to place ISIL on positive axis
            corrs(jj,mm)=-corrs(jj,mm);
        end
        ev1(jj,mm)=evdiag(1,1);
        [~,assorts(jj,mm)]=assortscalar(simnet,x);
        simvec=symadj2vec(simnet);
        lsqerr(jj,mm)=sqrt(mean((simvec-obsvec)'*(simvec-obsvec)));
    end
    nsigCorr=2*min(length(find(corrs(:,mm)<rObs)),...
        length(find(corrs(:,mm)>=rObs)));
    pvalCorr(mm)=nsigCorr/nrun;
    nsigAssort=2*min(length(find(assorts(:,mm)<thetaObs)),...
        length(find(assorts(:,mm)>=thetaObs)));
    pvalAssort(mm)=nsigAssort/nrun;
    rSimAv(mm)=mean(corrs(:,mm));
    rSimDev(mm)=std(corrs(:,mm));
    thetaSimAv(mm)=mean(assorts(:,mm));
    thetaSimDev(mm)=std(assorts(:,mm));
    ev1Av(mm)=mean(ev1(:,mm));
    netMedian(:,:,mm)=median(netall(:,:,:,mm),3);
end

if M==1
    summTab=table(scalName,rObs,pvalStand,rSimAv,rSimDev,pvalCorr,...
        thetaObs,thetaSimAv,thetaSimDev,pvalAssort,ev1Av);
else
    summTab=[];
end


lsqerrAv=mean(lsqerr);
corrAv=mean(corrs);
assortAv=mean(assorts);

[errBest,iBest]=min(lsqerrAv);
netBest=round(netav(:,:,iBest));
lamBest=lam(iBest);
assortBest=assortAv(iBest);
corrBest=corrAv(iBest);

[~,~,U,S]=modularity_corr(netBest,x,eigOrder,identeigvec,modmatrix);

outputFile=[outputFile '_sum'];
save(outputFile,'pars','corrs','assorts','pvalCorr','pvalAssort',...
    'rObs','thetaObs','ev1','netav','pvalStand','nodelabels','x',...
    'lsqerr','lsqerrAv','summTab','adjmat','corrAv','assortAv',...
    'lamBest','iBest','netBest','errBest','lamBest','assortBest',...
    'corrBest')

if ~isempty(summTab)
    writetable(summTab,[outputFile '.xlsx'])
end


[~,~,U,S]=modularity_corr(netBest,x,eigOrder,identeigvec,modmatrix);

if U(ixISIL,1)<0
    U(:,1)=-U(:,1);  % place ISIL on right
end

d1=U(:,1);
d2=U(:,2);
d3=U(:,3);


axsign=[1 1];
d1=axsign(1)*d1;
d2=axsign(2)*d2;

if identeigvec==1
    identdim=d1;
    compdim=d2;
    xnetlabel='Eigenvector 1';
    ynetlabel='Eigenvector 2';
elseif identeigvec==2
    identdim=d2;
    compdim=d1;
    xnetlabel='Eigenvector 2';
    ynetlabel='Eigenvector 1';
elseif identeigvec==3
    identdim=d3;
    compdim=d1;
    xnetlabel='Eigenvector 3';
    ynetlabel='Eigenvector 1';
end




idTypeFull=idNamesFull{idIdx};
idType=idNames{idIdx};
[cmax,pval]=corr(identdim,x);


% calculate assortativity
[~,assortcoef]=assortscalar(netBest,x);

nGroups=length(nodelabels);


% calculate classification p-val via binomial dist
if ~isempty(find(strcmpi('Mid',idNames),1))
    ixJihad=x>=mean(x);
    
else
    ixJihad=x>=3;
end
classPos=identdim>=0;
classNeg=identdim<0;
numPosRight=length(find(ixJihad==classPos));
numNegRight=length(find(ixJihad==classNeg));

if numPosRight>numNegRight
    classBest=classPos;
    numRight=numPosRight;
else
    classBest=classNeg;
    numRight=numNegRight;
end
pvalClass=binocdf(numRight,nGroups,0.5,'upper');
[pfit,pci]=binofit(numRight,nGroups);
missed=nodelabels(ixJihad~=classBest);
missedText='';
for kk=1:length(missed)
    if kk~=length(missed)
        missedText=[missedText missed{kk} ','];
    else
        missedText=[missedText missed{kk}];
    end
end

nodeColor=repmat([0 0 0],nGroups,1);
nodeColor(ixJihad,:)=repmat([1 1 1],length(find(ixJihad)),1);


%     nodeColor=repmat([0 0 1],nGroups,1);
axlabs={xnetlabel,ynetlabel};
figtitle=idTypeFull;
if binNetForPlot==1
    binmat=netBest;
    binmat(netBest>0)=1;
end
net2dplot(identdim,compdim,binmat,nodeColor,nodelabels,figtitle,axlabs)
ax=axis;



text(0.2*(ax(2)-ax(1))+ax(1),0.4*(ax(4)-ax(3))+ax(3),...
    {['assort = ' num2str(assortcoef,'%8.4f')]
    ['r = ' num2str(cmax,'%8.4f')]
    ['pst = ' num2str(pval,'%8.4f')]
    ['\lambda_1 = ' num2str(S(1,1),'%8.4f')]
    ['\lambda_2 = ' num2str(S(2,2),'%8.4f')]
    ['pcl = ' num2str(pvalClass,'%8.4f')]
    ['\gamma = ' num2str(pfit,'%6.3f') ' (' num2str(pci(1),'%6.3f')...
    ',' num2str(pci(2),'%6.3f') ')']
    ['missed: ' missedText]
    ['error: ' num2str(lsqerrAv(iBest),'%8.4f')]
    },'fontsize',10);
hold on
line([0 0],[ax(3) ax(4)],'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
hold off
netfig=gcf;
labelfig(netfig,[outputFile '_bestnet'])
saveas(netfig,[outputFile '_bestnet'])

% errfig=figure;
% plot(lam,lsqerrAv,'-o')
% xlabel('lam');ylabel('mean error')
% title(idTypeFull)
% labelfig(errfig,[outputFile '_err_' idType])
% saveas(errfig,[outputFile '_err_' idType])
% 
% corrfig=figure;
% plot(lam,corrAv,'-o')
% xlabel('lam');ylabel('mean correlation')
% title(idTypeFull)
% labelfig(corrfig,[outputFile '_corr_' idType])
% saveas(corrfig,[outputFile '_corr_' idType])
% 
% assortfig=figure;
% plot(lam,assortAv,'-o')
% xlabel('lam');ylabel('mean assortativity')
% title(idTypeFull)
% labelfig(assortfig,[outputFile '_assort_' idType])
% saveas(assortfig,[outputFile '_assort_' idType])

combofig=figure;
subplot(3,1,1)
plot(lam,corrAv)
hold on
plot([lam(1) lam(end)],[abs(corrObs) abs(corrObs)],'--')
hold off
currAx=gca;
currAx.XTickLabel=[];
ax=axis;
text(0.7*(ax(2)-ax(1))+ax(1),0.8*(ax(4)-ax(3))+ax(3),...
    {['corr = ' num2str(corrBest,'%8.4f')]})
ylabel('mean corr.')
title(idTypeFull)

subplot(3,1,2)
plot(lam,assortAv)
hold on
plot([lam(1) lam(end)],[assortObs assortObs],'--')
hold off
currAx=gca;
currAx.XTickLabel=[];
ax=axis;
text(0.7*(ax(2)-ax(1))+ax(1),0.8*(ax(4)-ax(3))+ax(3),...
    {['assort = ' num2str(assortBest,'%8.4f')]})
ylabel('mean assort.')

subplot(3,1,3)
plot(lam,lsqerrAv)
ax=axis;
text(0.7*(ax(2)-ax(1))+ax(1),0.8*(ax(4)-ax(3))+ax(3),...
    {['lam = ' num2str(lamBest,'%8.4f')]})
xlabel('lam');ylabel('mean error')
labelfig(combofig,[outputFile '_combo'])
saveas(combofig,[outputFile '_combo'])
