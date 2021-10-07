%%%Cost function for the curve fitting program: minimises the cost over
%%%steady state and time courses

function y=CostFunIAAUpGFPDownAlpha2(par)

%Load data
load('IAAUpGFPDownAlphaFull3')
% load('ParameterEvalIAAUpGFP')
TS=1000;
%Run simulation to compute cost function
y0=0;
par=par.^2;

for j=1:length(Input)
    [~,yODE1]=ode15s(@ActiveIAAUpGFPDownAlpha,[0 TS],zeros(5,1),odeset('refine',10),0,0,par); %Simulating the system for different input
    [t,yODE]=ode15s(@ActiveIAAUpGFPDownAlpha,[Time(1) TS],yODE1(end,:),odeset('refine',10),Input(j,1),Input(j,2),par);
    yEnd=yODE(:,end);        
    Y=binlin(t,yEnd,Time);
    K=3;
    K1=0; %Weight parameter for steady state values
    K2=0; %Weight parameter for second point (speed of response)
    y0=y0+sum((log(Y)-log(Data(:,j))).^2)+K1*((log(yEnd(end))-log(Data(end,j)))^2)+K2*(log(Y(2))-log(Data(2,j)))^2; %Cost function contribution of the first time series    
    if j==8 || j==13 || j==10  || j==11% || j==5
        y0=y0+K*sum((log(Y)-log(Data(:,j))).^2)+K1*((log(yEnd(end))-log(Data(end,j)))^2)+K2*(log(Y(2))-log(Data(2,j)))^2; %Cost function contribution of the first time series
    end
    if j==13 % || j==5
        y0=y0+sum((log(Y)-log(Data(:,j))).^2)+K1*((log(yEnd(end))-log(Data(end,j)))^2)+K2*(log(Y(2))-log(Data(2,j)))^2; %Cost function contribution of the first time series
    end
%     if j==3
%         AA=Y;
%     elseif j==9
%         BB=Y;
%     end
    %y0=y0+sum(abs(log(Y)-log(Data(:,j))));%+abs(log(yEnd(end))-log(Data(end,j)));
end    

% y0=y0+2*sum((log(AA)-log(Y)).^2);

y=sqrt(y0);

end