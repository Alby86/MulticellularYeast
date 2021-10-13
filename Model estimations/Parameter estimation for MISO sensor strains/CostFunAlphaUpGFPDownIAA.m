%%%Cost function for the curve fitting program: minimises the cost over
%%%steady state and time courses

function y=CostFunAlphaUpGFPDownIAA(par)

%Load data
load('AlphaUpGFPDownIAA')
TS=1000;
%Run simulation to compute cost function
y0=0;
[~,yODE1]=ode15s(@ActiveAlphaUpGFPDownIAA,[0 TS],zeros(5,1),odeset('refine',10),0,0,par); %Simulating the system for the initial state
    
for j=1:length(Input)
    [t,yODE]=ode15s(@ActiveAlphaUpGFPDownIAA,[Time(1) TS],yODE1(end,:),odeset('refine',10),Input(j,1),Input(j,2),par);
    yEnd=yODE(:,end);        
    Y=binlin(t,yEnd,Time);
    y0=y0+sum((log(Y)-log(Data(:,j))).^2); %Cost function contribution of the first time series    
    %y0=y0+sum(abs(log(Y)-log(Data(:,j))));%+abs(log(yEnd(end))-log(Data(end,j)));
end

y=sqrt(y0);

end