function [OptimalRatio,Y]=OptimizeGate(StrainCombo,Sensor,logic,N)

%Define the optimization parameters
A=zeros(2*(N-1),N);
for i=1:N-1 %to compute the bounds on the stochiometry parameters
    A(2*(i-1)+1,1)=1;
    A(2*i,1)=-1;
    A(2*i-1,i+1)=-10;
    A(2*i,i+1)=0.1;
end
b=zeros(2*(N-1),1);
Aeq=[];
beq=[];
lb=0.1*ones(N,1);
ub=10*ones(N,1);
options=optimoptions('fmincon','TolFun',1e-9,'TolCon',1e-9,'MaxIter',10000);
K0=ones(1,N); %Initial stochiometry

%Define Simulation time, parameters and inputs
load('ParameterListF')
Tt=1000; %time for transients
TS=12; %time for effective simulation
IAAMax=1000000;
alphaMax=100000;
betaMax=10000;

%Computing the background for each sensor
par0=ones(20,1);
Input0=zeros(1,10);
Output0=zeros(1,10);
input0=[0 0 0];
inputM=[alphaMax IAAMax betaMax];
x0=zeros(28,1);
Trans=1;
[~,yODE2]=ode15s(@TwoStrainSimulatorRatio,[0 Tt],x0,odeset('refine',10),input0,par0,Input0,Output0,K0,Trans);
X0=yODE2(end,:);
X0(3)=0; 
X0(6)=0;
X0(X0<0)=0;
[~,y1]=ode15s(@TwoStrainSimulatorRatio,[0 TS],X0,odeset('refine',10),input0,par0,Input0,Output0,K0,0);
[~,y2]=ode15s(@TwoStrainSimulatorRatio,[0 TS],X0,odeset('refine',10),inputM,par0,Input0,Output0,K0,0);
[~,y3]=ode15s(@TwoStrainSimulatorRatio,[0 TS],X0,odeset('refine',10),[0 IAAMax 0],par0,Input0,Output0,K0,0);
[~,y4]=ode15s(@TwoStrainSimulatorRatio,[0 TS],X0,odeset('refine',10),[alphaMax 0 0],par0,Input0,Output0,K0,0);
Back=[y2(end,9) y2(end,12) y1(end,15) y1(end,18) y4(end,23) y3(end,28)];

%Loading parameters
load('ParameterListF')

%Running Optimization
for i=1:size(StrainCombo,1)    
    
    %Assigning strain inputs and outputs
    parTot=[];
    inputTot=[];
    outputTot=[];
    for j=1:N
        parTot=[parTot par(StrainCombo(i,j),:)];
        inputTot=[inputTot input(StrainCombo(i,j),:)];
        outputTot=[outputTot output(StrainCombo(i,j),:)];
    end
    
    [Cost0(i),Y(:,2*(i-1)+1)]=Logic(K0,parTot,inputTot,outputTot,logic(StrainCombo(i,end),:),StrainCombo(i,end-1),Sensor,Back(Sensor),N);
    LogicHelper = @(x)Logic(x,parTot,inputTot,outputTot,logic(StrainCombo(i,end),:),StrainCombo(i,end-1),Sensor,Back(Sensor),N);
    [K1(i,:),Cost1(i)]=fmincon(LogicHelper,K0,A,b,Aeq,beq,lb,ub,[],options);
    [~,Y(:,2*i)]=Logic(K1(i,:),parTot,inputTot,outputTot,logic(StrainCombo(i,end),:),StrainCombo(i,end-1),Sensor,Back(Sensor),N);

end

OptimalRatio=[K1 Cost0' Cost1'];

end