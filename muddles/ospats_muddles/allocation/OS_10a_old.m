function [strat0,StratBest,ObarBest,ObarH1,d2,d20] = OS_10a(x,y,z_pred,s2,R2,ra,H,ncy,start,maxrun,d20)
% SCRIPT FOR OPTIMAL STRATIFICATION BY ITERATIVE RE-ALLOCATION. VERSION 6-11-2013
% Algorithm main features: (1) minimisation by 'first opportinity decent', 
% not 'steepest decent', (2) outer loop through grid points, not strata, 
% (3) calculation of reductions of objective function by re-allocations via
% sums of squared distances, not repeated function evaluation.
% Version 2 Sept 2013 - BM Error correction
% Version 3 Oct 2013 - JdG Multiple runs and 2 different initial solutions
% Version 5 Nov 2013 - JdG Prediction error variance spatially variable
% Version 8 5/12/2013 - JdG Different maximum covariance per pair of points

%                            INPUT
% x,y,z,s2: coordinates, prediction, and variance of prediction
% si, ra: sill and range corresponding to variance
% H : number of strata
% ncy : maximum number of cycles per run
% maxrun : number of runs

%if(nargin < 6), ncy = 100; end;   
n = length(x);

% generate initial solution
if start == 1  % use Cum-sqrt-f
    nclass = 100;
    [xm, strat0] = cumsqf(z, nclass, H);
else  % use kmeans 
    strat0 = kmeans([x,y],H,'replicates',10,'MaxIter',300);
end %if start==1
    
%                            INITIATION
% Vectorised calculation of distance matrix, without loop

if(d20==0),
% calculate distance matrix
xy = [x,y];
x1 = permute(xy, [ 1 3 2 ]);
y1 = permute(xy, [ 3 1 2 ]);
Lag = sqrt(sum((x1(:,ones(1,n),:) - y1(ones(1,n),:,:)).^2, 3));

% calculate Z difference (vectorized)
z1 = permute(z_pred, [ 1 3 2 ]);
z2 = permute(z_pred, [ 3 1 2 ]);
Z2 = (sum((z1(:,ones(1,n),:) - z2(ones(1,n),:,:)).^2, 3));
Z2 = Z2/R2;         %compensation for leveling

% calculate variance 
zs1 = permute(s2, [ 1 3 2 ]);
zs2 = permute(s2, [ 3 1 2 ]);
S2 = sum((zs1(:,ones(1,n),:) + zs2(ones(1,n),:,:)), 3); % S2i + S2j

% calculate matrix of maximum of covariance
    Cov_max = 0.5*S2;
   Cov = Cov_max.*exp(-3*Lag./ra);
   d2 = Z2 + S2 -2*Cov;
   %d2 = S2 -2*Cov; % emulate equal z-pred's
   %d2 = Lag;    %use only Lag
   %d2 = Z2;     % use only Z (true or pred)
    d20=d2;
else
    d2=d20;
end

format compact
%MinZ2 = min(min(Z2)), MeanZ2 = mean(mean(Z2)), MaxZ2 = max(max(Z2)); 
%MinS2 = min(min(S2)), MeanS2 = mean(mean(S2)), MaxS2 = max(max(S2));
%MinCov = min(min(Cov)), MeanCov = mean(mean(Cov)), MaxCov = max(max(Cov));

TOTd2 = sum(sum(d2))/2;
ObjNoStr = sqrt(TOTd2);
ObarH1 = ObjNoStr/n

% Calculate sums of d2's within strata (Sd2) and contributions of strata to Obj (cbObj).
Sd2 = zeros(1,H);
for strat = 1:H
    Sd2(strat) = 0;
    for i=1:(n-1)
        if strat0(i) == strat
            for j=(i+1):n
                if strat0(j) == strat                
                    Sd2(strat) = Sd2(strat)+d2(i,j);
                end %if strat0(j) == strat 
            end %for j=(i+1):n
        end %if strat0(i) == strat
    end %for i=1:(n-1)   
end %for strat = 1:H

Sd2_init = Sd2;

cbObj = sqrt(Sd2);
Obj = sum(cbObj);
ObarInit =Obj/n

%                           START RUNS
ObarBest = 1000;        %arbitrary large number to start with
StratBest = zeros(n);

for run = 1:maxrun
    TotTransf=0;
    Sd2 = Sd2_init;
    disp(['-----------------------------------RUN NR. ', num2str(run)])

    %disp('--------------start transferring-----------”')
    %                             ITERATIVE RE-ALLOCATION
    stratcy = strat0;       % stores stratum number for each point
    change = 1;

for cycle = 1:ncy        % loop through cycles
    transfers = 0;
    u = randperm(n);	% put grid points in random order

for t = u           % loop through grid points
    Delta = 0;
	change = 0;           % indicator of transfer
    A = stratcy(t);
    % remove t from A
    ij = find(stratcy == A);
    dA = sum(d2(t,ij))-d2(t,t);
    sumd2tinA = dA;   
    Sd2Amint = Sd2(A)-sumd2tinA;
    cbObjA = sqrt(Sd2Amint);
        
    for stratnr = 1:H   % idem between t and the points in different strata
        Delta = 0;
        sumd2plus = 0;
        if stratnr ~= A 
            % Add to B
            B = stratnr;
            ij = find(stratcy == B);
            dB = sum(d2(t,ij));
            sumd2plus = dB;
           
            cbObjB = sqrt(Sd2(B) + sumd2plus);
            Delta = cbObjA + cbObjB - cbObj(A) - cbObj(B);
            if Delta < -Obj*1e-10           % realize transfer
                change = 1;
                transfers = transfers+1;
                stratcy(t) = B;    % update stratification 
                Sd2(A) = Sd2(A)-sumd2tinA;
                Sd2(B) = Sd2(B)+sumd2plus;
                cbObj = sqrt(Sd2);
                Obj = sum(cbObj);
                Delta = 0;
            end %if Delta<0
            %if change == 1, bre}ak, end
        end %if stratnr~=A 
        if change == 1, break, end
    end %for stratnr=1:H  
end %for t=u

disp(['Cycle nr.: ', num2str(cycle), '   transfers:  ', num2str(transfers)])

TotTransf = TotTransf + transfers;
if (transfers == 0), break, end; % stopping rule
end %cy=1:ncy

%disp('--------------end transferring--------')
disp(['Number of cycles= ', num2str(cycle)])
disp(['Total number of transfers= ', num2str(TotTransf)])
Obj = sum(cbObj);
ObarFinal = Obj/n

% Update if last ObarFinal is smaller than the smallest previous ObarFinal
if ObarFinal < ObarBest
    ObarBest = ObarFinal
    StratBest = stratcy;
end %if ObarFinal < ObarBest

end %for run=1:maxrun

% final output
%Initial_distr = tabulate(strat0)
%Final_distr = tabulate(StratBest)













