function main()

    % ==========
    % Parameters
    % ==========

    % =====
    % Stock
    % =====

    sigma=.3;	% Volatility

    Su=2;    	% Upper Bound in price grid
    Sb=.00001;  % Lower Bound in price grid. Not = 0 to avoid S=0 in Black-Scholes Formula

    % =============
    % Interest rate and Dividend Yield
    % =============

    r=.03;
    D=.04;

    % =======
    % Options
    % =======

    T=1;	% Time to Maturity
    K=1;	% Strike Price
    KC=.8;
    KP=1.2;

    % ===========
    % Grid set up
    % ===========

    n=100; % number of points on the grid (including first and last)
    I=35;  % number of points on the grid (including first and last)


    % ================
    % Creation of Grid
    % ================

    dS=(Su-Sb)/(I-1); % Price step
    dt=T/(n-1);			% Time step

    v1=dt/dS^2;			% Two quantities needed in the Finite Differece Equation
    v2=dt/dS;


    SS=[Sb:dS:Su]';  	% Price Intervals
    TT=[0:dt:T]; 		% Time Intervals

    % =============================================
    % Initialization of V(S,k) and other quantities
    % =============================================

    %V=zeros(I,n);				% By definition, they are all matrices with
    Delta=zeros(I,n);			% dimension equal to the one of the grid
    Theta=zeros(I,n);
    Gamma=zeros(I,n);
    %VP=zeros(I,n);
    VC=zeros(I,n);

    BB=NaN*ones(I,n);
    EP=zeros(I,n);
    % =================================
    % Final conditions (over i for k=1)
    % =================================

    disp('	Computing V(S,T) using Final Conditions')

    for i=1:I
       %V(i,1)=max(SS(i)-KC,0)+max(KP-SS(i),0);

       %VP(i,1)=max(KP-SS(i),0);
       VC(i,1)=max(SS(i)-KC,0);
    end


    % ===================
    % Boundary Conditions 
    % ===================

        % Over k for i=1
    disp('	Computing V(0,t) using Lower Boundary Condition')

    for k=2:n
        %V(1,k)=KP;   %*exp(-r*(k-1)*dt);
        %VP(1,k)=KP;
        VC(1,k)=0;
    end

       % Over k for i=I

    disp('	Computing V(Su,t) using Upper Boundary Condition')

    for k=2:n
        %V(I,k)=SS(I)-K*exp(-r*(k-1)*dt);    	  % Notice: the (k-1) because since we cannot 
                                                             % have k=0  we must have as first element 
                                                              % k=1 (which corresponds to zero time 
                                                                  % to maturity

        %V(I,k)=SS(I)-KC;
        %VP(I,k)=0;
        VC(I,k)=SS(I)-KC;
    end


    % =====================   
    % Loop for other points
    % =====================

    disp('	Start Looping to Compute V(i,k): Please wait ')

    for k=2:n;  % Loop over time to maturity
        for i=2:I-1 % Loop over i

            a=1/2*sigma^2*SS(i)^2;   % Black-Scholes a(S,t)
            %a=1/2*(sigma-.3769+exp(-SS(i)))^2*SS(i)^2;  % Changing VOlatility Case
            b=(r-D)*SS(i);
            c=-r;

            A=v1*a-1/2*v2*b;
            B=-2*v1*a+dt*c;
            C=v1*a+1/2*v2*b;

            %ER(i,k)=A*VP(i-1,k-1)+(1+B)*VP(i,k-1)+C*VP(i+1,k-1);
            %V(i,k)=max(max(A*V(i-1,k-1)+(1+B)*V(i,k-1)+C*V(i+1,k-1),KP-SS(i)),SS(i)-KC);         
            %VP(i,k)=max(A*VP(i-1,k-1)+(1+B)*VP(i,k-1)+C*VP(i+1,k-1),KP-SS(i));
            VC(i,k)=max(A*VC(i-1,k-1)+(1+B)*VC(i,k-1)+C*VC(i+1,k-1),SS(i)-KC);

            %if VP(i,k)==KP-SS(i)
            %   BB(i,k)=1;
            %end
        end
    end


end
