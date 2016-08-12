function [klas,Obar] = OS_10a_all(x,y,z,s2,R2,ra,H, StratBest, d2,Sd2, xt,yt,zt,st2)
% allocation to existing OSPATS class
%
%                            INPUT
% x,y,z,s2: coordinates, prediction, and variance of prediction for input
% grid
% si,ra: sill and range corresponding to variance
% H : number of strata
% Stratbest: strata of inputs
% d2: distance ^2 matrix 
% Sd2: sum of distance matrix
% xt,yt,zt,st2: datum for allocation
% Output:
%   klas= allocated class
%   Obar  = calculated Obj function
n = length(x);

% distance of new point t to existing data
    Zt2=(z-zt).^2;
    St2=s2+st2;
    Lagt=sqrt((x-xt).^2+(y-yt).^2);          
    Zt2 = Zt2/R2;         %compensation for leveling

% calculate matrix of maximum of covariance
   Cov_max = 0.5*St2;
   Cov = Cov_max.*exp(-3.*Lagt./ra);
   dt2 = Zt2 + St2 -2.*Cov;

   clear Zt2 St2 Lagt Cov x y z s2;

% Insert the new point to d2
    d2(n+1,n+1)=0;
    d2(n+1,1:n)=dt2';
    d2(1:n,n+1)=dt2;
    A=1;    % try with strata 1
    n=n+1;
    StratBest(n)= A;
 % Calculate sums of d2's within strata (Sd2) and contributions of strata to Obj (cbObj).
     ij = find(StratBest == A);
     dA = sum(d2(n,ij));
     Sd2(A)=Sd2(A)+dA;   
     cbObj = sqrt(Sd2);

for strat = 2:H

    % remove t from A
    A=strat-1;
    t=n;
    ij = find(StratBest == A);
    dA = sum(d2(t,ij))-d2(t,t);
    sumd2tinA = dA;   
    Sd2Amint = Sd2(A)-sumd2tinA;
    cbObjA = sqrt(Sd2Amint);

    % add to B
    B = strat;
    ij = find(StratBest == B);
    dB = sum(d2(t,ij));
    sumd2plus = dB;
           
    cbObjB = sqrt(Sd2(B) + sumd2plus);
    Delta = cbObjA + cbObjB - cbObj(A) - cbObj(B);
 
    if Delta < 0           % Onbj function decrease? then  transfer
        StratBest(t) = B;    % update stratification 
        Sd2(A) = Sd2(A)-sumd2tinA;
        Sd2(B) = Sd2(B)+sumd2plus;
        cbObj = sqrt(Sd2);
        Obj = sum(cbObj);
        Delta = 0;
    end %if Delta<0
            
end %for strat = 1:H

klas=StratBest(t);
Obj = sum(cbObj);
Obar = Obj/n;













