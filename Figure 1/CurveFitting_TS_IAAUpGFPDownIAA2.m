clear all
% close all
clc

%Code to estimate ODE parameters
load('IAAUpGFPDownAlphaFull2')

load('ParameterEvalIAAUpGFPDownAlpha8');
parff=parf;
% parff=sqrt(parf);
clear parf parF Err

%First: Hill function steady state parameters
nopar=16; %Number of parameters

%Define the constrains for the optimization problem
A=-eye(nopar);
A(1,2)=1;
A(3,4)=1;
% A(5,7)=1;
b=zeros(nopar,1);
Aeq = [];
beq = [];
lb=zeros(nopar,1);
ub=Inf*ones(nopar,1);
ub(6)=5; %All the parameters are bound above by infinity, but the Hill coefficient, which is bounded by 5
ub(10)=5;
% lb(6)=1;
% lb(8)=1;
%ub(1:4)=1e-4;
% ub(nopar)=min(min(Data))*0.9;
% ub(10:11)=1;

%Then we start the optimization program
Rep=10; %Number of different initial conditions for parameter estimation
Err0=10000*sum(sum(Data));

for i=1:Rep
    i
    try
%     if mod(i,2)==0
        par=parff.*(1+rand(size(parff)).*sign(rand(size(parff))-0.5));
%     else
%         par(1:4)=10^(-4-round(rand*5))*rand(4,1);%Initial guess for the parameters
%         par(5:nopar)=rand(nopar-4,1);
%         par(6)=ub(6)*rand;
%         par(8)=ub(8)*rand;
%     end
%     par=parff.*(1+rand(size(parff)));
%     par(5:nopar)=rand(nopar-4,1);
%     par=rand(nopar,1);
    %Optimization and cost function
    options=optimoptions('fmincon','TolFun',1e-8,'TolCon',1e-8,'MaxIter',3000);
%     [parF(i,:),Err(i)]=fmincon(@CostFunIAAUpGFPDownAlpha,par,A,b,Aeq,beq,lb,ub,[],options);
    [parF(i,:),Err(i)]=fminsearch(@CostFunIAAUpGFPDownAlpha2,sqrt(par));
    if Err(i)<Err0 %If we found a new minimum of the cost function
        Err0=Err(i);
        parf=parF(i,:);
    end      
    catch
    end
end

parf=parf.^2;

TS=1000;
figure
[~,yODE1]=ode15s(@ActiveIAAUpGFPDownAlpha,[0 TS],zeros(5,1),odeset('refine',10),0,0,parf); %Simulating the system for different input
for j=1:length(Input)
    [t,yODE]=ode15s(@ActiveIAAUpGFPDownAlpha,[Time(1) Time(end)],yODE1(end,:),odeset('refine',10),Input(j,1),Input(j,2),parf);
    yEnd=yODE(:,end);
    %Y=binlin(t,yEnd,Time);
    plot(Time,Data(:,j),'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
    hold on
    plot(t,yEnd,'LineWidth',2,'DisplayName','Simulation')
end 


filename='ParameterEvalIAAUpGFPDownAlpha9';
save(filename,'parf','parF','Err')