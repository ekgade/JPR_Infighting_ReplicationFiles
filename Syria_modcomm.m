function Syria_modcomm(pars)

outFileFull=pars.outFileFull;
eigOrder=pars.eigOrder;
numnets=pars.numnets;
binarize=pars.binarize;
binNetForPlot=pars.binNetForPlot;
identeigvec=pars.identeigvec;
numeigs=3;

adjmat=pars.net;
nodelabels=pars.netgroups;
if numnets==2
    adjmat2=pars.netred2;
    nodelabels2=pars.netgroups;
end
nettype=pars.nettype;

idNames=pars.idNames;
idNamesFull=pars.idNamesFull;
idScores=pars.idScores;


if binarize
    adjmat(adjmat>0)=1;
    outFileFull=[outFileFull '_bin'];
else
    outFileFull=[outFileFull '_wei'];
end
degs=sum(adjmat,2);

izdeg=find(degs==0);
if ~isempty(izdeg)
    degs(izdeg)=[];
    nodelabels(izdeg)=[];
    adjmat(izdeg,:)=[];
    adjmat(:,izdeg)=[];
    idScores(izdeg,:)=[];
end


[eigVec,eigVal]=modularityeigvec(adjmat,numeigs,eigOrder);

outFileFull=[outFileFull '_ev' num2str(identeigvec)];


if identeigvec==1
    identdim=eigVec(:,1);
    compdim=eigVec(:,2);
    xnetlabel='Eigenvector 1';
    ynetlabel='Eigenvector 2';
elseif identeigvec==2
    identdim=eigVec(:,2);
    compdim=eigVec(:,1);
    xnetlabel='Eigenvector 2';
    ynetlabel='Eigenvector 1';
elseif identeigvec==3
    identdim=eigVec(:,3);
    compdim=eigVec(:,2);
    xnetlabel='Eigenvector 3';
    ynetlabel='Eigenvector 1';
end

ixISIL=find(strcmp('ISIL',nodelabels),1); % place ISIL on positive side
if identdim(ixISIL)<0
    identdim=-identdim;
end



numId=length(idNames);
for idIndex=1:numId
    tp=idScores(:,idIndex);
    idType=idNames{idIndex};
    idTypeFull=idNamesFull{idIndex};
    [cmax,pval]=corr(identdim,tp);
%     tpfit = linspace(min(tp),max(tp),100);
%     bz = polyfit(tp,identdim,1);
%     sdfit = polyval(bz,tpfit);

    idfit = linspace(min(identdim),max(identdim),100);
    bz = polyfit(identdim,tp,1);
    sdfit = polyval(bz,idfit);

    % calculate assortativity
    [~,assortcoef]=assortscalar(adjmat,tp);
    
    nGroups=length(nodelabels);
    
    
    % calculate classification p-val via binomial dist
    ixJihad=tp>=3;
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


[~,degassort]=assortscalar(adjmat,degs);


    
    
    axlabs={xnetlabel,ynetlabel};
    figtitle=[nettype ': ' idType];
    if binNetForPlot==1
        binadjmat=adjmat;
        binadjmat(adjmat>0)=1;
    end
    net2dplot(identdim,compdim,binadjmat,nodeColor,nodelabels,figtitle,axlabs)
    ax=axis;
    
        text(0.6*(ax(2)-ax(1))+ax(1),0.4*(ax(4)-ax(3))+ax(3),...
        {['assort = ' num2str(assortcoef,'%8.4f')]
        ['r = ' num2str(cmax,'%8.4f')]
        ['p = ' num2str(pval,'%8.4f')]
        ['\lambda_1 = ' num2str(eigVal(1),'%8.4f')]
        ['\lambda_2 = ' num2str(eigVal(2),'%8.4f')]
%         ['pcl = ' num2str(pvalClass,'%8.4f')]
%         ['\gamma = ' num2str(pfit,'%6.3f') ' (' num2str(pci(1),'%6.3f')...
%         ',' num2str(pci(2),'%6.3f') ')']
%         ['missed: ' missedText]
        },'fontsize',10);
    hold on
    line([0 0],[ax(3) ax(4)],'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
    hold off
    netfig=gcf;
    labelfig(netfig,[outFileFull '_net_' idType])
    saveas(netfig,[outFileFull '_net_' idType])

    axlabs={idTypeFull,xnetlabel};
    net2dplot(tp,identdim,[],nodeColor,nodelabels,figtitle,axlabs);
    hold on
    line(sdfit,idfit,'linestyle','-','linewidth',2,'color',[0.5 0.5 0.5])
    ax=axis;
    text(0.5*(ax(2)-ax(1))+ax(1),0.9*(ax(4)-ax(3))+ax(3),...
        {['r = ' num2str(cmax,'%8.4f') '; p = ' num2str(pval,'%8.4f')]
%         ['\gamma = ' num2str(pfit,'%6.3f') ' (' num2str(pci(1),'%6.3f')...
%         ',' num2str(pci(2),'%6.3f') ')']
%         ['missed: ' missedText]
        },'fontsize',10);
        idmid=3;

    line([idmid idmid],[ax(3) ax(4)],'color',[0.5 0.5 0.5],'linestyle','--',...
        'linewidth',2)
    line([ax(1) ax(2)],[0 0],'color',[0.5 0.5 0.5],'linestyle','--',...
        'linewidth',2)
    hold off
    corrfig=gcf;
    labelfig(corrfig,[outFileFull '_corr_' idType])
    saveas(corrfig,[outFileFull '_corr_' idType])
end

save(outFileFull,'pars','adjmat','nodelabels','idScores','eigVec',...
    'eigVal','degs','degassort')
    
%%%%% Second network %%%%%%

if numnets==2
    eigOrder2=pars.eigOrder2;
    identeigvec=pars.identeigvec2;
if binarize
    adjmat2(adjmat2>0)=1;
else
end
degs2=sum(adjmat2,2);

izdeg=find(degs2==0);
if ~isempty(izdeg)
    degs2(izdeg)=[];
    nodelabels2(izdeg)=[];
    adjmat2(izdeg,:)=[];
    adjmat2(:,izdeg)=[];
    idScores(izdeg,:)=[];
end


[eigVec2,eigVal2]=modularityeigvec(adjmat2,numeigs,eigOrder2);


if identeigvec==1
    identdim=eigVec2(:,1);
    compdim=eigVec2(:,2);
    xnetlabel='Eigenvector 1';
    ynetlabel='Eigenvector 2';
elseif identeigvec==2
    identdim=eigVec2(:,2);
    compdim=eigVec2(:,1);
    xnetlabel='Eigenvector 2';
    ynetlabel='Eigenvector 1';
elseif identeigvec==3
    identdim=eigVec2(:,3);
    compdim=eigVec2(:,2);
    xnetlabel='Eigenvector 3';
    ynetlabel='Eigenvector 1';
end

[~,degassort2]=assortscalar(adjmat2,degs2);

numNodes=length(nodelabels2);
nodeColor=repmat([0 0 0],numNodes,1);
    
    
    axlabs={xnetlabel,ynetlabel};
    figtitle=pars.nettype2;
    idType2=[idType '2'];
    if binNetForPlot==1
        binadjmat=adjmat2;
        binadjmat(adjmat2>0)=1;
    end
    net2dplot(identdim,compdim,binadjmat,nodeColor,nodelabels2,figtitle,axlabs)
    ax=axis;
    

    
    text(0.2*(ax(2)-ax(1))+ax(1),0.4*(ax(4)-ax(3))+ax(3),...
        {['deg. assort = ' num2str(degassort2,'%8.4f')]
        ['\lambda_1 = ' num2str(eigVal2(1),'%8.4f')]
        ['\lambda_2 = ' num2str(eigVal2(2),'%8.4f')]
        },'fontsize',10);
    hold on
    line([0 0],[ax(3) ax(4)],'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5]);
    hold off
    netfig=gcf;
    labelfig(netfig,[outFileFull '_net_' idType2])
    saveas(netfig,[outFileFull '_net_' idType2])
end
    

