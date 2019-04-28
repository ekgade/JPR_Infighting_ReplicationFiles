function pars=Syria_modcomm_pars(tPeriod,iVar,identEV,testType)

% tPeriod: data time period ('Full','2014','2015)
% iVar: independent variable ('ideo','power')
% identEV: test eigenvector (1,2)
% testType: homophily or heterophily ('hom','het')


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

pars.identeigvec=identEV;

if strcmp(testType,'hom')
    pars.eigOrder='la';
elseif strcmp(testType,'het')
    pars.eigOrder='sa';
end



netDataFile2='';


delCore={};   % do not delete any groups
% delCore={'ISIL'};  % uncomment to delete ISIL



rootFolder='.';

outFolder=[rootFolder filesep 'results'];



outFileFull=[outFolder filesep outFileBase];
DataFolder=[rootFolder filesep 'input'];



pars.binarize=0;
pars.modmatrix='modularity';

pars.binNetForPlot=1;

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

% if isequal(pars.nettype,'joint ops')
%     pars.eigOrder='la';
% elseif isequal(pars.nettype,'infighting')
%     pars.eigOrder='sa';
% end


if pars.numnets==2
    if isequal(pars.nettype2,'joint ops')
    pars.eigOrder2='la';
elseif isequal(pars.nettype2,'infighting')
    pars.eigOrder2='sa';
    end
end


if pars.numnets==2
    pars.identeigvec2=1;
end

idDataFileFull=[DataFolder filesep idDataFile];
load(idDataFileFull)

pars.tgrps=tgrps;
pars.idScoresOrig=idScores;
pars.idNames=idNames;
pars.idNamesFull=idNamesFull;

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





