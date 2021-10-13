clear all
close all
clc

%Code to estimate ODE parameters for MISO strain
load('AlphaUpGFPDownIAA')

nopar=16; %Number of parameters

%Define the constrains for the optimization problem
A=-eye(nopar);
b=zeros(nopar,1);
Aeq = [];
beq = [];
lb=zeros(nopar,1); %All parameters are forced to be bigger or equal to 0
ub=Inf*ones(nopar,1);
ub(6)=5; %All the parameters are bound above by infinity, but the Hill coefficient, which is bounded by 5
ub(10)=5; 

%Then we start the optimization program
Rep=10; %Number of different initial conditions for parameter estimation
Err0=10000*sum(sum(Data));

for i=1:Rep

    par=rand(nopar,1);%Initial guess for the parameters
    
    %Optimization and cost function
    options=optimoptions('fmincon','TolFun',1e-8,'TolCon',1e-8,'MaxIter',3000);
    [parF(cont,:),Err(cont)]=fmincon(@CostFunAlphaUpGFPDownIAA,par,A,b,Aeq,beq,lb,ub,[],options);
    if Err(i)<Err0 %If we find a new minimum of the cost function
        Err0=Err(cont);
        parf=parF(cont,:);
    end      
end

TS=1000;
figure
[~,yODE1]=ode15s(@ActiveAlphaUpGFPDownIAA,[0 TS],zeros(5,1),odeset('refine',10),0,0,parf); %Simulating the system for different input
for j=1:length(Input)
    [t,yODE]=ode15s(@ActiveAlphaUpGFPDownIAA,[Time(1) Time(end)],yODE1(end,:),odeset('refine',10),Input(j,1),Input(j,2),parf);
    yEnd=yODE(:,end);
    %Y=binlin(t,yEnd,Time);
    plot(Time,Data(:,j),'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
    hold on
    plot(t,yEnd,'LineWidth',2,'DisplayName','Simulation')
end 

set(gca, 'XScale','log')
set(gca,'FontSize',15)
savefig('AlphaUpGFPDownIAATS')

filename='ParameterEvalAlphaUpGFPDownIAA';
save(filename,'parf','parF','Err')