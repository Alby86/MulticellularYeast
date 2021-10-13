%%%Cost function for the curve fitting program: minimises the cost over
%%%steady state and time courses

function y=CostFunAlphaUpIAA(par)

%Load data
load('AlphaUpIAA')
TS=1000;
par=par.^2;
K=1; %Ratio between the two strains
%Run simulation to compute cost function for the steady state values
[~,yODE1]=ode15s(@ActiveAlphaUpIAA,[0 TS],zeros(6,1),odeset('refine',10),0,par,K); %Simulating the system for different input
load('ParameterEvalIAADown2xGFP')
[~,yODE2]=ode15s(@Repress,[0 TS],zeros(3,1),odeset('refine',10),0,pa0rf);

for j=1:length(Input)
    [t,yODE]=ode15s(@ActiveAlphaUpIAA,[Time(1) Time(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),Input(j),par,K);
    yEnd(j)=yODE(end,end);
end
    
K1=0;
ySS=sum((log(yEnd)-log(Data')+K1*(log(yEnd(1))-log(Data(1)))).^2); %Cost function contribution of the first time series

%Now compute the Time Series simulation 
[t,yODE]=ode15s(@ActiveAlphaUpIAA,[TStime(1) TStime(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),TSinput,par);

yTS=binlin(t,yODE(:,end),TStime);

yTS=(sum(log(yTS)-log(TSdata))).^2;

y=sqrt(ySS+yTS);

end