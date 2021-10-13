clear all
close all
clc

%Code to estimate ODE parameters

%Load data
load('AlphaDownGFP')

nopar=9; %Number of parameters

%Define the constrains for the optimization problem
A=-eye(nopar);
b=zeros(nopar,1);
Aeq = [];
beq = [];
lb=zeros(nopar,1); %All parameters are forced to be bigger or equal to 0
ub=Inf*ones(nopar,1);
ub(5)=5; %All the parameters are bound above by infinity, but the Hill coefficient, which is bounded by 5

%Then we start the optimization program
%First, the first repeat
Rep=10; %Number of different initial conditions for parameter estimation
Err0=1e2*sum(sum(DataRep1));

for i=1:Rep
    i
    par=rand(nopar,1);%Initial guess for the parameters
    par(5)=rand;
    
    %Optimization and cost function
    options=optimoptions('fmincon','TolFun',1e-8,'TolCon',1e-8,'MaxIter',3000);
    [parF,Err]=fmincon(@CostFunAlphaDownGfpRep1,par,A,b,Aeq,beq,lb,ub,[],options);
    if Err<Err0 %If we found a new minimum of the cost function
        Err0=Err;
        parf1=parF;
    end      
    
end

%Then, the second repeat
Err0=1e2*sum(sum(DataRep2));

for i=1:Rep
    i
    par=parf1.*(ones(size(parf1))+rand(size(parf1)));
    %Optimization and cost function
    options=optimoptions('fmincon','TolFun',1e-8,'TolCon',1e-8,'MaxIter',3000);
    [parF,Err]=fmincon(@CostFunAlphaDownGfpRep2,par,A,b,Aeq,beq,lb,ub,[],options);
    if Err<Err0 %If we found a new minimum of the cost function
        Err0=Err;
        parf2=parF;
    end      
    
end

%Thirdly, the third repeat
Err0=1e2*sum(sum(DataRep3));

for i=1:Rep
    i
    par=parf2.*(ones(size(parf2))+rand(size(parf2)));
    %Optimization and cost function
    options=optimoptions('fmincon','TolFun',1e-8,'TolCon',1e-8,'MaxIter',3000);
    [parF,Err]=fmincon(@CostFunAlphaDownGfpRep3,par,A,b,Aeq,beq,lb,ub,[],options);
    if Err<Err0 %If we found a new minimum of the cost function
        Err0=Err;
        parf3=parF;
    end      
    
end

%Computing the average and standard deviation of the 3 parameter sets
parfAV=(parf1+parf2+parf3)./3;
parfSTD=std([parf1,parf2,parf3]');

%Now, computing the fitting for the average values
Err0=1e2*sum(sum(Average));

for i=1:Rep
    i
    par=parfAV.*(ones(size(parfAV))+rand(size(parfAV)));
    %Optimization and cost function
    options=optimoptions('fmincon','TolFun',1e-8,'TolCon',1e-8,'MaxIter',3000);
    [parF,Err]=fmincon(@CostFunAlphaDownGfp,par,A,b,Aeq,beq,lb,ub,[],options);
    if Err<Err0 %If we found a new minimum of the cost function
        Err0=Err;
        parf=parF;
    end      
    
end


%Plotting the average simulation vs the data
%Time Series
figure
TS=1000;
for j=1:length(Input)
    [~,yODE1]=ode15s(@Repress,[0 TS],[0 0 0],odeset('refine',10),0,parf); %Simulating the system for different input
    [t,yODE]=ode15s(@Repress,[Time(1) Time(end)],yODE1(end,:),odeset('refine',10),Input(j),parf); %Simulating the system for different input    
    y(:,j)=binlin(t,yODE(:,3),Time); %Resampling to the same time point as the data
    errorbar(Time,Average(:,j),Error(:,j),'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
    hold on
    plot(Time,y(:,j),'LineWidth',2,'DisplayName','Simulation')
end

set(gca,'FontSize',15)
savefig('AlphaDownGFPTS')

%Steady States
figure
TS=1000;
input=(Input(1):0.1:Input(end));
for j=1:length(input)
    [~,yODE1]=ode15s(@Repress,[0 TS],[0 0 0],odeset('refine',10),0,parf); %Simulating the system for different input
    [t,yODE]=ode15s(@Repress,[Time(1) Time(end)],yODE1(end,:),odeset('refine',10),input(j),parf); %Simulating the system for different input    
    Y(j)=yODE(end,3);
end
Input(1)=0.1; %To allow for logplots
errorbar(Input,Average(end,:),Error(end,:),'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
hold on
plot(input,Y,'LineWidth',2,'DisplayName','Simulation')

set(gca, 'XScale','log')
set(gca,'FontSize',15)

savefig('AlphaDownGFPSS')

%Saving the estimated parameter set
filename='ParameterEvalAlphaDownGFP';
save(filename,'parf','parf1','parf2','parf3','parfSTD','parfAV')