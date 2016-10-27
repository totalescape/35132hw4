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

sigma=0.73;	% Volatility
S0=7.85;

% =============
% Interest rate and Dividend Yield
% =============

r=.0;
D=0.0;

% =======
% Options
% =======

T=5;	% Time to Maturity
KC=20;	% Strike Price


% ===========
% Simulations
% ===========



dt = 1/(252*2);%1/252;
N=T/dt;
NSim=100000;

% Implementation
% ==============
dBt=sqrt(dt)*randn(NSim,N);
St=zeros(NSim,N);
St(:,1)=S0*ones(NSim,1);

% Conversion px
ConPx = NaN*zeros(NSim,N);
for t=2:N
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
MM(:,N)=max(SSit(:,N)-KC,0);
figure
for tt=N:-1:3; 
   disp('Time to Maturity')
   disp(1-tt/N)
   
   % Step 1: Select the path in the money at time tt-1
   
   I=find(SSit(:,tt-1)-KC >0 | SSit(:,tt)-KC >0);
   ISize=length(I);
   
   % Step 3: Project CashFlow at time tt onto basis function at time tt-1
   
   if tt==N
      %YY - Produces the PV of next period
      YY=(ones(ISize,1)*exp(-r*[1:N-tt+1]*dt)).*MM(I,tt:N);
      %YYFull=(ones(NSim,1)*exp(-r*[1:N-tt+1]*dt)).*MM(:,tt:N);
   else
      %YY - Produces the PV of all future periods (should be only 1
      %non-zero)
      YY=sum(((ones(ISize,1).*exp(-r*[1:N-tt+1]*dt)).*MM(I,tt:N))')';
      %YYFull=sum(((ones(NSim,1)*exp(-r*[1:N-tt+1]*dt)).*MM(:,tt:N))')';
   end

   %Perform regression using in the money data points
   SSb=SSit(I,tt-1);
   XX=[ones(ISize,1),SSb,SSb.^2,SSb.^3,SSb.^4,SSb.^5];
   BB=inv(XX'*XX)*XX'*YY;
   
   %Build regression data for later
   SSb2=SSit(:,tt-1);
   XX2=[ones(NSim,1),SSb2,SSb2.^2,SSb2.^3,SSb2.^4,SSb2.^5];
   
   plot(SSb,XX*BB,'.',SSb,SSb-KC,':'); %,SSb,YY,'*') %plot of (if exercise now value)
   %plot(SSb2,XX2*BB,'.',SSb2,SSb2-KC,':',SSb2,YYFull,'*') %plot of all
   legend('Expected Payoff if Wait','Payoff Today if Exercise')
   xlabel('Stock Price')
   title('Estimation of Exercise Frontier')
   pause(.0001)
   
   %                                 XX2*BB = Continuation value
   IStop=find(SSit(:,tt-1)-KC >= max(XX2*BB,0)); %IStop = Index of exercise now
   ICon=setdiff([1:NSim],IStop);               %ICon = Index of no-exercise
   
   MM(IStop,tt-1)=SSit(IStop,tt-1)-KC;         %Use exercise value for exercise now
   MM(IStop,tt:N)=zeros(length(IStop),N-tt+1); %Zero out previous exercises
   MM(ICon,tt-1)=zeros(length(ICon),1);        %Zero out non-exercise, keeping future exercise
   
   ErrorI = find(sum(MM(:,tt-1:N)')'<max(SSit(:,N)-KC,0));
   if (sum(ErrorI)>0) 
       temp = max(XX2*BB,0);
        [sum(MM(ErrorI,tt-1:N)')' max(SSit(ErrorI,N)-KC,0) SSit(ErrorI,tt-1)-KC temp(ErrorI) ErrorI SSit(ErrorI,tt-1) SSit(ErrorI,tt)]
   end
 end
 
 YY=sum(((ones(NSim,1)*exp(-r*[1:N-1]*dt)).*MM(:,2:N))')';
 
 EuroValue = mean(exp(-r*T).*max(SSit(:,N)-KC,0))

 Value=mean(YY);
 sterr=std(YY)/sqrt(NSim);
 
 disp('Value of American Put Option')
 disp(Value)
 disp('St. Error')
 disp(sterr)
 disp('Black')
 disp(Bsc(S0,KC,r,D,sigma,T))