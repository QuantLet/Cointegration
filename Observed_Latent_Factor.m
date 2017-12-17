
%    Stock and Watson (1988) "Testing for common trends" Journal of
%    American Statistical Association. Vol.83, 1097-1107
 
 
clear all;
clc;
tic
data=xlsread('CombinedResultsCDSIndices_IncludingDeterminantVariables.xls','CDSIndices2','B:AC');
 
level=data(:,[1 3 5 7 9 11 13 15]);
IG=data(:,[1 3 9 11]);
HY=data(:,[5 7 13 15]);
SR=data(:,[1 5 9 13]);
LR=data(:,[3 7 11 15]);
US=data(:,[1 3 5 7]);
EU=data(:,[9 11 13 15]);
 
%       tmpData=level;
%         tmpData=IG;
%        tmpData=HY;
%      tmpData=SR;
%      tmpData=LR;
%     tmpData=US;
%      tmpData=EU;
 
% PRe-Crisis
 t=130;  % setting in the table 5 results
 
% t=190;
 
   tmpData=level(1:t,1:8);
%         tmpData=IG(1:t,1:4);
%       tmpData=HY(1:t,1:4);
%     tmpData=SR(1:t,1:4);
%     tmpData=LR(1:t,1:4);
%     tmpData=US(1:t,1:4);
%     tmpData=EU(1:t,1:4);
 
% Post-Crisis (since Lehman Brother bankruptcy in Sep 2008) 
%    tmpData=level(t+1:283,1:8);    
%         tmpData=IG(t+1:283,1:4);
%       tmpData=HY(t+1:283,1:4);
%     tmpData=SR(t+1:283,1:4);
%     tmpData=LR(t+1:283,1:4);
%     tmpData=US(t+1:283,1:4);
%     tmpData=EU(t+1:283,1:4);
 
 
% change=data(:,[2 4 6 8 10 12 14 16]);
% IGchange=data(:,[2 4  10 12]);  % IG group
% HYchange=data(:,[6 8 14 16]);  % HY group
% ShortYchange=data(:,[2 6  10 14]);  % 5Y group
% LongYchange=data(:,[4 8 12 16]);  % 10Y group
% USchange=data(:,[2 4 6 8]);
% EUchange=data(:,[10 12 14 16]);
%      tmpData=change;
%      tmpData=IGchange;
%     tmpData=HYchange;
%   tmpData=ShortYchange;
%   tmpData=LongYchange;
%   tmpData=USchange;
%  tmpData=EUchange;
 
vstoxx=data(:,17);
vix=data(:,18);
aaa=data(:,19);
baa=data(:,20);
euro10Y=data(:,21);
euro1Y=data(:,22);
tbond10Y=data(:,23);
tbond5Y=data(:,24);
tbond1Y=data(:,25);
Libor=data(:,26);
Repo=data(:,27);
Tbill=data(:,28);
 
yieldcurve5.US=tbond5Y-tbond1Y;
yieldcurve10.US=tbond10Y-tbond1Y;
yieldcurve.EU=euro10Y-euro1Y;
CreditSpread=baa-aaa;
Counterparty=Libor-Repo;
Liquidity=Repo-Tbill;
 
 
dyieldcurve5.US(2:283)=(yieldcurve5.US(2:283)-yieldcurve5.US(1:282))*100;
dyieldcurve10.US(2:283)=(yieldcurve10.US(2:283)-yieldcurve10.US(1:282))*100;
 
dyieldcurve.EU(2:283)=(yieldcurve.EU(2:283)-yieldcurve.EU(1:282))*100;
dtbond1Y(2:283)=(tbond1Y(2:283)-tbond1Y(1:282))*100;
deuro1Y(2:283)=(euro1Y(2:283)-euro1Y(1:282))*100;
dvix(2:283)=(vix(2:283)-vix(1:282));
dvstoxx(2:283)=(vstoxx(2:283)-vstoxx(1:282));
dCreditSpread(2:283)=(CreditSpread(2:283)-CreditSpread(1:282))*100;
 
 
%  observed variables during pre-crisis
Y=[tbond1Y(1:t) CreditSpread(1:t) yieldcurve5.US(1:t) vix(1:t) Counterparty(1:t) Liquidity(1:t)];
 
%  observed variables during crisis
% Y=[tbond1Y(t+1:283) CreditSpread(t+1:283) yieldcurve5.US(t+1:283) vix(t+1:283) Counterparty(t+1:283) Liquidity(t+1:283)];
[nr,nl]=size(Y);
Ystat=zeros(nl,1);
Ypval=zeros(nl,1);
 
% for i=1:7
% [Ystat(i),Ypval(i)]=augdf(Y(:,i),0,4);
% end
 
 
[T,N]=size(tmpData);
[V,D]=eig(tmpData'*tmpData);
 
 
 k=5;    % number of factors for entire sample
%   k=3;      % number of factors for subsample
 
factor=tmpData*V(:,1:N)*sqrt(N)./N;
 
b=zeros(N,nl);
residual=zeros(T,nl);
Rstat=zeros(nl,1);
Rpval=zeros(nl,1);
R2=zeros(nl,1);
 
for i=1:nl
b(:,i)=inv(factor'*factor)*factor'*Y(:,i);
Yhat=factor*b(:,i);
% ADF test for residuals, if residual is I(0), then there exist
% cointegration even if F is I(1).
residual(:,i)=Y(:,i)-Yhat;
[Rstat(i),Rpval(i)]=augdf(residual(:,i),0,4);
R2(i)=var(Yhat)/var(Y(:,i));
 
 
end
