clear all
close all
clc

TS=1000;
%Code to estimate ODE parameters
load('BetaUpBAR1')

%Simulation of the initial conditions for the sensor strain
load('ParameterEvalAlphaUpGFP')
[~,yODE2]=ode15s(@Active,[0 TS],zeros(3,1),odeset('refine',10),0,parf);
clear parf

%Then, optimization problem definition
nopar=3; %Number of parameters

%Then we start the optimization program
Rep=10; %Number of different initial conditions for parameter estimation
Err0=1000*sum(sum(Data));

for i=1:Rep
    
    par=rand(nopar,1);%Initial guess for the parameters
    
    %Optimization and cost function
    options=optimoptions('fmincon','TolFun',1e-8,'TolCon',1e-8,'MaxIter',3000);
    [parF,Err]=fminsearch(@CostFunBetaUpBAR1,par);
    if Err<Err0 %If we found a new minimum of the cost function
        Err0=Err;
        parf=parF;
    end      
    
end

parf=parf.^2;

figure
TS=1000;

[~,yODE1]=ode15s(@ActivateBetaUpBAR1,[0 TS],zeros(6,1),odeset('refine',10),0,0,parf); %Simulating the system for different input
for i=1:length(Input1)
    Input2(1)=0;
    input2=(Input2(1):1:Input2(end));
    for j=1:length(input2)
        [t,yODE]=ode15s(@ActivateBetaUpBAR1,[Time(1) Time(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),input2(j),Input1(i),parf); %Simulating the system for different input    
        yEnd(j)=yODE(end,end);
    end
    Input2(1)=1;
    input2(1)=[];
    yEnd(1)=[];
    errorbar(Input2,Data(:,i),Error(:,i),'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
    hold on
    plot(input2,yEnd,'LineWidth',2,'DisplayName','Simulation')
    clear yEnd
end

set(gca, 'XScale','log')
set(gca,'FontSize',15)

savefig('BetaUpBAR1SS')

[t,yODE]=ode15s(@ActivateBetaUpBAR1,[Time(1) Time(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),TSinput(2),TSinput(1),parf); %Simulating the system for different input    

yTS=yODE(:,end);

figure
plot(TStime,TSdata,'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
hold on
plot(t,yTS,'LineWidth',2,'DisplayName','Simulation')

set(gca,'FontSize',15)

savefig('BetaUpBAR1TS')

filename='ParameterEvalBetaUpBAR1';
save(filename,'parf')