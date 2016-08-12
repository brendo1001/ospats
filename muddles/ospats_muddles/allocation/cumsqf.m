function [xm, strat]= cumsqf(z, nclass, H)
% cumsqf cumulative square root of f(y) stratification
% based on Delenius 
% input: x: data, nclass : number class in cumsqf calculation, ns : number
% of strata
% output xm, for each strata: (1) number of samples; (2) mean; (3) lower
% boundary (min); (4) upper boundary (max)

% calculate the distribution based on nclass
[d, io]= sort(z);
nd=size(d);
dmin=min(d);
dmax=max(d);
step=(dmax-dmin)/nclass;
h=[dmin:step:dmax];

for i=1:nclass,
    ll=h(i);
    ul=h(i+1);
    ij=find(d>=ll & d<ul);
    if(i==nclass), ij=find(d>=ll & d<=ul); end;
    nij=length(ij);
    fx(i,1)=ll;
    fx(i,2)=ul;
    fy(i,1)=nij;
end
sqf=sqrt(fy);
cumsqf=cumsum(sqf);

% set out boundaries based on ns (number of strata)
cmax=max(cumsqf);
cmin=min(cumsqf);
cstep=(cmax-cmin)./H;
cbound=[cmin:cstep:cmax];
% find out corresponding class
for i=2:H,
    cm=abs(cumsqf-cbound(i));
    [ym, ij]=min(cm);
    db(i-1)=ij;
end

idat=zeros(nd);
% calculate the strata on the actual data
for i=1:H,
    if(i==1), ll=1; else ll=db(i-1); end;
    if(i==H), ul=nclass; else ul=db(i); end;
    
    fll=fx(ll,1);
    ful=fx(ul,1);
    if(i==H), ful=fx(ul,2); end;
    
    ij=find(d>=fll & d<ful);
    if(i==H), ij=find(d>=fll & d<=ful); end;
    idat(ij)=i;
    dij=d(ij);
    nij=length(dij);
    
    xm(i,1)= nij;
    xm(i,2)=mean(dij);
    xm(i,3)=var(dij);
    xm(i,4)=min(dij);
    xm(i,5)=max(dij);
    
end

strat=ones(nd);
strat(io)=idat;
