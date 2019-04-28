function pars=Syria_netsim_pars(tPeriod,iVar,simType,nrun)


if strcmp(tPeriod,'Full')
    netDataFile='infight_full_Sep2017';
elseif strcmp(tPeriod,'2014')
    netDataFile='infight_2014_Sep2017';
elseif strcmp(tPeriod,'2015')
    netDataFile='infight_2015_Sep2017';
end

if strcmp(iVar,'ideo')
    outPrefix='Id';
    idDataFile='id_Sep2017';
elseif strcmp(iVar,'power')
    outPrefix='Pow';
    idDataFile='powmid_Sep2017';
end
outFileBase=[outPrefix tPeriod];

pars.simtype=simType;

if pars.simtype==0
    outFileBase=[outFileBase '_Null'];
elseif pars.simtype==1
    outFileBase=[outFileBase '_Like'];    
elseif pars.simtype==-1
    outFileBase=[outFileBase '_Anti'];
end

rootFolder='.';

outFolder=[rootFolder filesep 'results'];



outFileFull=[outFolder filesep outFileBase];
DataFolder=[rootFolder filesep 'input'];


delCore={};
% delCore={'ISIL'};  % uncomment to delete ISIL

if pars.simtype==0
    pars.lam=[];
else
    if strcmp(iVar,'ideo')
        pars.lam=[0.2:0.1:6]';
    elseif strcmp(iVar,'power')
        pars.lam=[500:500:25000]';
    end
end
pars.nrun=nrun;
pars.idIdx=1;
pars.binarize=0;
pars.modmatrix='modularity';
pars.binNetForPlot=1;
netDataFile2='';

netDataFileFull=[DataFolder filesep netDataFile];
load(netDataFileFull)
if exist('allgroups','var')
    pars.nodenames=allgroups;
    pars.net=jofullsym;
    pars.nettype='joint ops';
elseif exist('infallgroups','var')
    pars.nodenames=infallgroups;
    pars.net=infightfullsym;
    pars.nettype='infighting';
end

outFileFull=[outFileFull '_' pars.nettype];

if pars.nrun>1
    outFileFull=[outFileFull '_nr' int2str(pars.nrun)];
end

if ~isempty(netDataFile2)
    netDataFileFull2=[DataFolder filesep netDataFile2];
    load(netDataFileFull2)
    if isequal(pars.nettype,'infighting')
        pars.nodenames2=allgroups;
        pars.net2=jofullsym;
        pars.nettype2='joint ops';
    elseif isequal(pars.nettype,'joint ops')
        pars.nodenames2=infallgroups;
        pars.net2=infightfullsym;
        pars.nettype2='infighting';
    end
    [pars.netred1,pars.netred2,pars.netgroups]=...
        commonnets(pars.net,pars.net2,pars.nodenames,pars.nodenames2);
    pars.numnets=2;
else
    pars.netgroups=pars.nodenames;
    pars.nodenames2='';
    pars.net2=[];
    pars.numnets=1;
end

if isequal(pars.nettype,'joint ops')
    pars.eigOrder='la';
elseif isequal(pars.nettype,'infighting')
    pars.eigOrder='sa';
end

if pars.simtype==1  % if homophily is being tested
    pars.eigOrder='la';
end

if pars.simtype==-1  % if heterophily is being tested
    pars.eigOrder='sa';
end

if pars.numnets==2
    if isequal(pars.nettype2,'joint ops')
    pars.eigOrder2='la';
elseif isequal(pars.nettype2,'infighting')
    pars.eigOrder2='sa';
    end
end


pars.identeigvec=1;  % use 1st eigenvector for identity dimension
if pars.numnets==2
    pars.identeigvec2=1;
end

idDataFileFull=[DataFolder filesep idDataFile];
load(idDataFileFull)

pars.tgrps=tgrps;
pars.idScoresOrig=idScores;
pars.idNames=idNames;
pars.idNamesFull=idNamesFull;

if ~isempty(pars.lam)
   outFileFull=[outFileFull '_' idNames{pars.idIdx}]; 
end

grpNames=tgrps;


if ~isempty(delCore)
    outFileFull=[outFileFull '_no'];
end
ixDelCore=[];
for k=1:length(delCore)
    idx=find(strcmp(delCore{k},grpNames),1);
    if ~isempty(idx)
        ixDelCore=[ixDelCore;idx];
        outFileFull=[outFileFull '_' delCore{k}];
    end
end
grpNames(ixDelCore)=[];
idScores(ixDelCore,:)=[];

netgroups=pars.netgroups;
if pars.numnets==1
net=pars.net;
elseif pars.numnets==2
    net=pars.netred1;
    net2=pars.netred2;
end
% if IF deleted, add its ties to ASIM
if find(strcmp('IF',delCore),1)
    ixIF=find(strcmp('IF',netgroups),1);
    if ~isempty(ixIF)
        ixASIM=find(strcmp('ASIM',netgroups),1);
        net(ixASIM,:)=net(ixASIM,:)+net(ixIF,:);
        net(:,ixASIM)=net(ixASIM,:)';
        if pars.numnets==2
            net2(ixASIM,:)=net2(ixASIM,:)+net2(ixIF,:);
            net2(:,ixASIM)=net2(ixASIM,:)';
        end
        outFileFull=[outFileFull '_IF2ASIM'];
    end
end

if pars.binarize
    net(net>0)=1;
    outFileFull=[outFileFull '_bin'];
else
    outFileFull=[outFileFull '_wei'];
end


N=length(grpNames);
ixKeep=[];
ixAbsent=[];
for j=1:N
    idx=find(strcmp(grpNames{j},netgroups),1);
    if ~isempty(idx)
        ixKeep=[ixKeep;idx];
    else
        ixAbsent=[ixAbsent;j];
    end
end
ixDel=1:length(netgroups);
ixDel(ixKeep)=[];
tnet=net;
tnet(:,ixDel)=[];
tnet(ixDel,:)=[];
[~,keepOrd]=sort(ixKeep);
net=zeros(length(ixKeep));
net(:,keepOrd)=tnet;
net(keepOrd,:)=net;
net=net-diag(diag(net));
if pars.numnets==2
    tnet2=net2;
    tnet2(:,ixDel)=[];
    tnet2(ixDel,:)=[];
    net2=zeros(length(ixKeep));
    net2(:,keepOrd)=tnet2;
    net2(keepOrd,:)=net2;
    net2=net2-diag(diag(net2));
end
pars.netfull=pars.net;
pars.net=net;
if pars.numnets==2
    pars.net2full=pars.net2;
    pars.net2=net2;
end

pars.ixAbsent=ixAbsent;
pars.grpsAbsent=grpNames(ixAbsent);
grpNames(ixAbsent)=[];
idScores(ixAbsent,:)=[];

pars.netgroups=grpNames;
pars.idScores=idScores;

pars.outFileFull=outFileFull;





