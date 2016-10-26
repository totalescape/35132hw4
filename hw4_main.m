% TN 11: American Options via Monte Carlo Simulations
% ===================================================
close all
clear all

%tic % starting the clock timer

% ==========
% Parameters
% ==========

% =====
% Stock
% =====

sigma=.3;	% Volatility
S0=1; 36;

% =============
% Interest rate and Dividend Yield
% =============

r=.03;
D=.04;

% =======
% Options
% =======

T=5;	% Time to Maturity
KP=1.2;	% Strike Price


% ===========
% Simulations
% ===========



dt =1/50;
N=T/dt;
NSim=20000;


% Implementation
% ==============
dBt=sqrt(dt)*randn(NSim,N);
St=zeros(NSim,N);
St(:,1)=S0*ones(NSim,1);

   for t=2:N;t
      St(:,t)=St(:,t-1).*exp((r-D-.5*sigma^2)*dt+sigma*dBt(:,t)); 
      
      
      
   end
   
   SSit=St;   
   

NSim=size(SSit,1);

% Computing the value of American Options
% =======================================

% Work Backwards

% Initialize CashFlow Matrxi
%N=N+1;

MM=NaN*ones(NSim,N);
MM(:,N)=max(KP-SSit(:,N),0);
figure
for tt=N:-1:3; 
   disp('Time to Maturity')
   disp(1-tt/N)
   
   % Step 1: Select the path in the money at time tt-1
   
   I=find(KP-SSit(:,tt-1)>0);
   ISize=length(I);
   
   % Step 3: Project CashFlow at time tt onto basis function at time tt-1
   
   if tt==N
   YY=(ones(ISize,1)*exp(-r*[1:N-tt+1]*dt)).*MM(I,tt:N);
	else
      YY=sum(((ones(ISize,1)*exp(-r*[1:N-tt+1]*dt)).*MM(I,tt:N))')';
   end

   SSb=SSit(I,tt-1);
   XX=[ones(ISize,1),SSb,SSb.^2,SSb.^3,SSb.^4,SSb.^5];
   BB=inv(XX'*XX)*XX'*YY;
   
   
   
   
   SSb2=SSit(:,tt-1);
   XX2=[ones(NSim,1),SSb2,SSb2.^2,SSb2.^3,SSb2.^4,SSb2.^5];
   
   plot(SSb,XX*BB,'.',SSb,KP-SSb,':')
   legend('Expected Payoff if Wait','Payoff Today if Exercise')
   xlabel('Stock Price')
   title('Estimation of Exercise Frontier')
   pause(.0001)
   
   
   IStop=find(KP-SSit(:,tt-1)>=max(XX2*BB,0));
   ICon=setdiff([1:NSim],IStop);
   
   MM(IStop,tt-1)=KP-SSit(IStop,tt-1);
   MM(IStop,tt:N)=zeros(length(IStop),N-tt+1);
   MM(ICon,tt-1)=zeros(length(ICon),1);
   
 end
 
  YY=sum(((ones(NSim,1)*exp(-r*[1:N-1]*dt)).*MM(:,2:N))')';

 Value=mean(YY);
 sterr=std(YY)/sqrt(NSim);
 
 disp('Value of American Put Option')
 disp(Value)
 disp('St. Error')
 disp(sterr)