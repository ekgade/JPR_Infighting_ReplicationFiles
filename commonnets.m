function [acom,bcom,labs]=commonnets(a,b,alab,blab)

[alabsort,ixasort]=sortrows(alab);
asort=a(ixasort,:);
asort=asort(:,ixasort);
[blabsort,ixbsort]=sortrows(blab);
bsort=b(ixbsort,:);
bsort=bsort(:,ixbsort);


labs=intersect(alabsort,blabsort);
[~,ixua,ixub]=setxor(alabsort,blabsort);

acom=asort;
acom(ixua,:)=[];
acom(:,ixua)=[];

bcom=bsort;
bcom(ixub,:)=[];
bcom(:,ixub)=[];

