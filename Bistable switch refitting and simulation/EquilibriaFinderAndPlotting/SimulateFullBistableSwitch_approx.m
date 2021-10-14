clear all
close all
clc

%Load final parameter values
K=[0.2 42.5 25 0.7]; %For K4=0.8, the 'high' state becomes unstable
SHIFTIAA=200; %Exogenous IAA added to the media
SHIFTalpha=0;

%Define the IAA and alpha pulses parameters
DilRate=3/2; %Rate of dilution (2 means that it is diluted of 1/2 
DilTime=40; %How often the dilution happens
T1=3*60; %Time of the exogenous IAA induction pulse
T2=15*60; %Time of the exogenous alpha induction pulse
U2=50000; %IAA peak max intensity
U1=50000; %alpha peak max intensity

%Simulation parameters
Ts=60*48;
Tt=100;
u=[000 0000]; %Inputs= [alpha IAA]

%Initial conditions
load('ParameterEvalIAADown2xGFP')
[~,y0]=ode15s(@Repress,[0 Tt],[0 0 0],odeset('refine',10),0,parf);  
x0=[0 0 y0(end,:)]; %Initial conditions [alpha0, IAA0, sensor(0)] 
%Option with the whole system in coculture
% [~,y0]=ode15s(@BistableSwitchFullTS,[0 Ts],zeros(1,5),odeset('refine',10),0,K,0,0);  
% x0=y0(end,:);

%Simulating the system
[t,yODE]=ode15s(@BistableSwitchFullTS,[0 Ts],x0,odeset('refine',10),u,K,SHIFTIAA,SHIFTalpha);  

%Now, plotting the sensor profile with background adjustment
load('BackGroundAdj','Correction')
figure
subplot(1,2,1)
plot(t/60,yODE(:,end)*Correction-1,'LineWidth',2,'DisplayName','Simulation')
set(gca,'FontSize',15)
xlabel('Time(hrs)');
ylabel('Relative Mean Fluorescence');
box('off')

%And the experimental data on top
hold on
load('BistableExpData')
errorbar(data4st(:,1),data4st(:,2)*Correction-1,data4st(:,3)*Correction,'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
legend('show')
xlim([0 32])

%The IAA pulse profile
subplot(1,2,2)
for i=1:length(t)
    if t(i)>T1
        Dil11(i)=floor((t(i)-T1)/DilTime);
        uts2(i)=U2/(DilRate^Dil11(i));
    else
        uts2(i)=0;
    end
end
plot(t/60,uts2,'LineWidth',2,'DisplayName','IAA pulse')
set(gca,'FontSize',15)
hold on

%The alpha pulse profile
for i=1:length(t)
    if t(i)>T2
        Dil12(i)=floor((t(i)-T2)/DilTime);
        uts1(i)=U1/(DilRate^Dil12(i));
    else
        uts1(i)=0;
    end
end
plot(t/60,uts1,'LineWidth',2,'DisplayName','Alpha pulse')
set(gca,'FontSize',15)
xlabel('Time(hrs)');
ylabel('Concentration (nM)');
box('off')
legend('show')

savefig('SymBistableSwitchApprox_split')

%Now a figure with the sensor profile and the pulsed on the same graph
figure
plot(t/60,yODE(:,end)*Correction-1,'LineWidth',2,'DisplayName','Simulation')
set(gca,'FontSize',15)
xlabel('Time(hrs)');
ylabel('Relative Mean Fluorescence');
box('off')

%And the experimental data on top
hold on
load('BistableExpData')
errorbar(data4st(:,1),data4st(:,2)*Correction-1,data4st(:,3)*Correction,'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
xlim([0 32])

%And the input peaks
plot(t/60,uts2/1e3,'LineWidth',2,'DisplayName','IAA pulse')
plot(t/60,uts1/1e3,'LineWidth',2,'DisplayName','Alpha pulse')

legend('show')

savefig('SymBistableSwitchApprox_All')

save('SymBistableSwitchApprox')