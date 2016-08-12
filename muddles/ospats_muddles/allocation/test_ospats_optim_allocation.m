%%%%%%%%%%%%%%
% OC surface
clear all
load Nowley2014.txt

x=Nowley2014(:,1);
y=Nowley2014(:,2);
z=Nowley2014(:,3);
se=Nowley2014(:,4);
s2=se.^2.*2;
space=10;   % grid spacing


% convert to grid index
xmin=min(x);
ymin=min(y);
xmax=max(x);
ymax=max(y);
nr=round((ymax-ymin+0.5*space)/space)+1;
nc=round((xmax-xmin+0.5*space)/space)+1;
yi=round((y-ymin+0.5*space)/space)+1;
xj=round((x-xmin+0.5*space)/space)+1;

% sample every 10th grid spacing
jj=find(mod(xj,10) ==1 & mod(yi,10) ==1);
xs=x(jj);
ys=y(jj);
zs=z(jj);
s2s=s2(jj);
d20=0;
scatter(xs,ys,[],zs, 'filled'); %%scatter plot
exDat= [xs,ys, zs, s2s];
%x=xs;
%y=ys;
%z=zs;
%s2=s2s;

% Run Ospats
ra=4000;
ncy=400;
start=2;
maxrun=1;
R2= 0.42;
H=5;
d20=0;
[strat0,StratBest,ObarBest,ObarH1,d2,Sd2] = OS_10a(xs,ys,zs,s2s,R2,ra,H,ncy,start,maxrun,d20);
figure(1);scatter(xs,ys,[],StratBest, 'filled');


% Allocate to a New grid
jk=find(mod(xj,3) ==1 & mod(yi,3) ==1);
%jk=[1:length(x)]';
ij=ismember(jk,jj); % find if there are data from Ostapts grid
im=find(ij==0);

xst=x(im);yst=y(im);zst=z(im);s2st=s2(im);
exDat1= [xst,yst, zst, s2st]
dlmwrite('myFile_all.txt',exDat1)
tic
strat_t=zeros(length(xst),1);
Ov=zeros(length(xst),1);
for it=1:length(xst) % allocate point by point
    it
    xt=xst(it);
    yt=yst(it);
    zt=zst(it);
    st2=s2st(it);
    [klas,Obar] = OS_10a_all(xs,ys,zs,s2s,R2,ra,H, StratBest, d2,Sd2, xt,yt,zt,st2);    
    strat_t(it)=klas;
    Ov(it)=Obar;
    klas
end
toc
figure(2);
scatter([xs;xst],[ys;yst],[],[StratBest;strat_t], 'filled');

xp=[1:length(xst)];np=length(xp);
plot(xp,Ov, xp,ones(np,1).*ObarBest);

