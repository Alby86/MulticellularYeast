%%%Cost function for the curve fitting program: minimises the cost over
%%%steady state and time courses

function y=CostFunBetaUpBAR1(par)

%Load data
load('BetaUpBAR1')
TS=1000;
par=par.^2;
%Run simulation to compute cost function
[~,yODE1]=ode15s(@ActivateBetaUpBAR1,[0 TS],zeros(6,1),odeset('refine',10),0,0,par); %Simulating the system for different input
load('ParameterEvalAlphaUpGFP')
[~,yODE2]=ode15s(@Active,[0 TS],zeros(3,1),odeset('refine',10),0,parf);

for i=1:length(Input1)
    for j=1:length(Input2)
        [t,yODE]=ode15s(@ActivateBetaUpBAR1,[Time(1) Time(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),Input2(j),Input1(i),par);
        yEnd(j,i)=yODE(end,end);
    end
end
    
ySS=sum(sum((log(yEnd)-log(Data)).^2)); %Cost function contribution of the first time series

[t,yODE]=ode15s(@ActivateBetaUpBAR1,[Time(1) Time(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),TSinput(2),TSinput(1),parf); %Simulating the system for different input   

yTS=binlin(t,yODE(:,end),TStime);

yTS=(sum(log(yTS)-log(TSdata))).^2;

y=sqrt(ySS+yTS);

end