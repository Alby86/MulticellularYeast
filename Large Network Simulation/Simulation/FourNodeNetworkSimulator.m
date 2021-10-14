%Code to simulate any 4 random strain combination from the library

clear all
close all
clc

%Simulation time and and inputs
Tt=1000; %time for transients
TS=12; %time for effective simulation
IAA=[0 100 500 1000 5000];
alpha=[0 1 5 10 50 200 1000];
beta=[0 1 5 10 100];
%for now, no IAM as input

Kratio=[0.1 0.5 1 5 10]; %Ratios between the two strains

load('ParameterListF')
N=size(par,1);
x0=zeros(34,1);
inputT=zeros(1,size(input,2)); %input array for the transient (no cell-to-cell communication)
outputT=zeros(1,size(output,2));
cont=0;
cont2=0;
cont3=0;
Trans=1;
Sim=0;
%Space pre-allocation
SpN=length(IAA)*length(alpha)*length(beta)*N*(N-1)*(N-2)*(N-3)/24;
SpNSS=length(IAA)*length(alpha)*N*(N-1)*(N-2)*(N-3)/24;
Y1=zeros(SpN,TS+1);
Y2=zeros(SpN,TS+1);
Y3=zeros(SpN,TS+1);
Y4=zeros(SpN,TS+1);
Y5=zeros(SpN,TS+1);
Y6=zeros(SpN,TS+1);
Y1SS=zeros(SpNSS,length(beta));
Y2SS=zeros(SpNSS,length(beta));
Y3SS=zeros(SpNSS,length(beta));
Y4SS=zeros(SpNSS,length(beta));
Y5SS=zeros(SpNSS,length(beta));
Y6SS=zeros(SpNSS,length(beta));
Yindex=zeros(SpN,7);
YSSindex=zeros(SpNSS,6);

%Build the simulation network
for i=1:N
    
    for j=i+1:N
        
        for kk=j+1:N
            
            for qq=kk+1:N
        
        cont3=cont3+1;
        input0=[0,0,0];
        parS=[par(i,:),par(j,:),par(kk,:),par(qq,:)];
        inputS0=[inputT,inputT,inputT,inputT];
        outputS0=[outputT,outputT,outputT,outputT];
        K0=ones(4,1);
        [~,yODE2]=ode15s(@FourStrainSimulatorRatio,[0 Tt],x0,odeset('refine',10),input0,parS,inputS0,outputS0,K0,Trans);
        X0=yODE2(end,:);
        X0(3)=0; 
        X0(6)=0;
        x0(31)=0;
        x0(34)=0;
        %for loops on the input concentration
        for k=1:length(IAA)
            for q=1:length(alpha)
                cont2=cont2+1;
                for w=1:length(beta)
                    cont=cont+1;
                    INPUTS=[alpha(q),IAA(k),beta(w)];
                    inputS=[input(i,:),input(j,:),input(kk,:),input(qq,:)];
                    outputS=[output(i,:),output(j,:),output(kk,:),output(qq,:)];
                    [t,yODE]=ode15s(@FourStrainSimulatorRatio,[0 TS],X0,odeset('refine',10),INPUTS,parS,inputS,outputS,K0,Sim);
                    Y1(cont,:)=binlin(t,yODE(:,9),(0:TS)); %alpha-factor sensor output
                    Y1SS(cont2,w)=yODE(end,9); %Steady State for the alpha-factor sensor                   
                    Y2(cont,:)=binlin(t,yODE(:,12),(0:TS)); %IAA sensor output
                    Y2SS(cont2,w)=yODE(end,12); %Steady State for the IAA sensor 
                    Y3(cont,:)=binlin(t,yODE(:,15),(0:TS)); %alpha-factor upregulated sensor output
                    Y3SS(cont2,w)=yODE(end,15); %Steady State for the alpha-factor upregulated sensor
                    Y4(cont,:)=binlin(t,yODE(:,18),(0:TS)); %IAA upregulated sensor output
                    Y4SS(cont2,w)=yODE(end,18); %Steady State for the IAA upregulated sensor
                    Y5(cont,:)=binlin(t,yODE(:,23),(0:TS)); %MISO IAA->GFP|-alpha sensor output
                    Y5SS(cont2,w)=yODE(end,23); %Steady State for the MISO IAA->GFP|-alpha sensor
                    Y6(cont,:)=binlin(t,yODE(:,28),(0:TS)); %MISO alpha->GFP|-IAA sensor output
                    Y6SS(cont2,w)=yODE(end,28); %Steady State for the MISO alpha->GFP|-IAA sensor
                    Yindex(cont,:)=[i,j,kk,qq,IAA(k),alpha(q),beta(w)]; %Index for the time series
                    YSSindex(cont2,:)=[i,j,kk,qq,IAA(k),alpha(q)];
                    Ylabel={'Strain #1','Strain #2','IAA (nM)','alpha (nM)','beta (nM)'}; 
                    
                    if k==1 && q==1 && w==1 %no input
                        GATElabel(cont3,1:4)=[i,j,kk,qq];
                        GATE(cont3,1:6)=[yODE(end,9) yODE(end,12) yODE(end,15) yODE(end,18) yODE(end,23) yODE(end,28)];
                    elseif k==1 && q==1 && w==5 %max beta
                        GATE(cont3,7:12)=[yODE(end,9) yODE(end,12) yODE(end,15) yODE(end,18) yODE(end,23) yODE(end,28)];
                    elseif k==1 && q==7 && w==1 %max alpha
                        GATE(cont3,13:18)=[yODE(end,9) yODE(end,12) yODE(end,15) yODE(end,18) yODE(end,23) yODE(end,28)];
                    elseif k==1 && q==7 && w==5 %max alpha&beta
                        GATE(cont3,19:24)=[yODE(end,9) yODE(end,12) yODE(end,15) yODE(end,18) yODE(end,23) yODE(end,28)];
                    elseif k==5 && q==1 && w==1 %max IAA
                        GATE(cont3,25:30)=[yODE(end,9) yODE(end,12) yODE(end,15) yODE(end,18) yODE(end,23) yODE(end,28)];
                    elseif k==5 && q==1 && w==5 %max IAA&beta
                        GATE(cont3,31:36)=[yODE(end,9) yODE(end,12) yODE(end,15) yODE(end,18) yODE(end,23) yODE(end,28)];
                    elseif k==5 && q==7 && w==1 %max IAA&alpha
                        GATE(cont3,37:42)=[yODE(end,9) yODE(end,12) yODE(end,15) yODE(end,18) yODE(end,23) yODE(end,28)];
                    elseif k==5 && q==7 && w==5 %max IAA&alpha&beta
                        GATE(cont3,43:48)=[yODE(end,9) yODE(end,12) yODE(end,15) yODE(end,18) yODE(end,23) yODE(end,28)];
                    end
                end
            end
        end
            end
        end
    end
end

%plotting the t-SNE separation
% y=tsne(Y1,Yindex(:,1:2),2,2,30);

% y=tsne(Y1,[],2,2,30);
% gscatter(y(:,1), y(:,2),Yindex(:,1:2));
  
filename='FourNodeNetworkSimulation';
save(filename)