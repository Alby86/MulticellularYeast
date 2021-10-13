clear all
close all
clc


%Code to estimate ODE parameters
load('AlphaUpIAA')

%Simulation of the initial conditions for the sensor strain
load('ParameterEvalIAADown2xGFP')
Tt=1000;
[~,yODE2]=ode15s(@Repress,[0 Tt],zeros(3,1),odeset('refine',10),0,parf);
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
    [parF,Err]=fminsearch(@CostFunAlphaUpIAA,par);
    if Err<Err0 %If we found a new minimum of the cost function
        Err0=Err;
        parf=parF;
    end      
    
end

parf=parf.^2; %The estimated parameters were squared during optimization to consider only positive parameters

figure
K=1; %Ratio between the two strains
Tt=1000;
M=max(max(Input,ValidationInput));
input=(Input(1):1:M);
[~,yODE1]=ode15s(@ActiveAlphaUpIAA,[0 Tt],zeros(6,1),odeset('refine',10),0,parf,K); %Simulating the system for different input
for j=1:length(input)
    [t,yODE]=ode15s(@ActiveAlphaUpIAA,[Time(1) Time(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),input(j),parf,K); %Simulating the system for different input    
    yEnd(1,j)=yODE(end,end);
end

K=5; %Ratio between the two strains
Tt=1000;
M=max(max(Input,ValidationInput));
input=(Input(1):1:M);
[~,yODE1]=ode15s(@ActiveAlphaUpIAA,[0 Tt],zeros(6,1),odeset('refine',10),0,parf,K); %Simulating the system for different input
for j=1:length(input)
    [t,yODE]=ode15s(@ActiveAlphaUpIAA,[Time(1) Time(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),input(j),parf,K); %Simulating the system for different input    
    yEnd(2,j)=yODE(end,end);
end

K=10; %Ratio between the two strains
Tt=1000;
M=max(max(Input,ValidationInput));
input=(Input(1):1:M);
[~,yODE1]=ode15s(@ActiveAlphaUpIAA,[0 Tt],zeros(6,1),odeset('refine',10),0,parf,K); %Simulating the system for different input
for j=1:length(input)
    [t,yODE]=ode15s(@ActiveAlphaUpIAA,[Time(1) Time(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),input(j),parf,K); %Simulating the system for different input    
    yEnd(3,j)=yODE(end,end);
end

Input(1)=1;
input(1)=[];
yEnd(:,1)=[];
errorbar(Input,Data,Error,'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
hold on
for j=1:3
    plot(input,yEnd(j,:),'LineWidth',2,'DisplayName','Simulation')
end
errorbar(ValidationInput,ValidationData,ValidationError,'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')

set(gca, 'XScale','log')
set(gca,'FontSize',15)

savefig('AlphaUpIAASS')

K=1;
[t,yODE]=ode15s(@ActiveAlphaUpIAA,[TStime(1) TStime(end)],[yODE1(end,1:2) 0 yODE2(end,1:3)],odeset('refine',10),TSinput,parf,K); %Simulating the system for different input    
yTS=yODE(:,end);

figure
plot(TStime,TSdata,'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
hold on
plot(t,yTS,'LineWidth',2,'DisplayName','Simulation')

set(gca,'FontSize',15)

savefig('AlphaUpIAATS')

filename='ParameterEvalAlphaUpIAA';
save(filename,'parf')