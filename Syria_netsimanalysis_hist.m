function Syria_netsimanalysis_hist(pars,R,L,alpha)

if nargin<4
    alpha=.05;
end
inputFile=[pars.outFileFull '_sum'];

idNamesFull=pars.idNamesFull;
idIdx=pars.idIdx;


load(inputFile,'lsqerr')
outputFile=[inputFile '_hist_' 'R_' int2str(R) '_L_' int2str(L)];

lam=pars.lam;
nrun=pars.nrun;
M=length(lam);

errMin=zeros(R,1);
lamMin=zeros(R,1);

a=zeros(L,M);

for j=1:R
    for k=1:M
        a(:,k)=randsample(lsqerr(:,k),L,true);
    end
    [errMin(j),idx]=min(mean(a));
    lamMin(j)=lam(idx);
end
  
errMinAv=mean(errMin);
lamMinAv=mean(lamMin);

ci=prctile(lamMin,[100*alpha/2,100*(1-alpha/2)]);

f1=figure;
hErr=histogram(lamMin);
ax=axis();
a1=gca;
hold on
plot(a1,[lamMinAv lamMinAv],[ax(3) ax(4)],'color','k','LineWidth',1);
plot(a1,[ci(1) ci(1)],[ax(3) ax(4)],'color','r','LineWidth',1);
plot(a1,[ci(2) ci(2)],[ax(3) ax(4)],'color','r','LineWidth',1);
textstats={
    ['mean=' num2str(lamMinAv)]
    ['ci=' num2str(ci)]
    };
xpos=0.65*(ax(2)-ax(1))+ax(1);
ypos=0.75*(ax(4)-ax(3))+ax(3);
text(xpos,ypos,textstats)
title(idNamesFull{idIdx})
hold off

labelfig(f1,outputFile)
saveas(f1,outputFile)



save(outputFile,'pars','errMin','errMinAv','lamMin','lamMinAv','ci','hErr')


