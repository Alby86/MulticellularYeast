%%%Cost function for the curve fitting program: minimises the cost over
%%%steady state and time courses

function y=BistableSwitchFullTS(t,x,w,K,SHIFTIAA,SHIFTalpha)

%Load parameters
load('ParSScon_1Arm_V5')
parf1=parf;
clear parf PAR ERR
load('ParSScon_2Arm_V5')
parf2=parf;
clear parf PAR ERR
load('ParSScon_2Arm_V5WithExIAA')
parff2=parf;
clear parf PAR ERR
load('ParSScon_1Arm_V5WithExIAA')
parff1=parf;
clear parf PAR ERR
load('Arm1TimePar','parf')
parf1TS=parf;
clear parf PAR ERR
load('Arm2TimePar','parf')
parf2TS=parf;
clear parf PAR ERR

%Reproducing the exogenous input conditions as in the experiment
DilRate=3/2; %Rate of dilution (2 means that it is diluted of 1/2) 
DilTime=40; %How often the dilution happens
T1=3*60; %Time of the exogenous IAA induction pulse
T2=15*60; %Time of the exogenous alpha induction pulse
U2=50000; %IAA peak max intensity
U1=75000; %alpha peak max intensity
if t>T1 %IAA induction
    Dil1=floor((t-T1)/DilTime);
    u(2)=U2/(DilRate^Dil1);
else
    u(2)=0;
end
if t>T2 %alpha induction
    Dil2=floor((t-T2)/DilTime);
    u(1)=U1/(DilRate^Dil2);
else
    u(1)=0;
end

%Defining ODEs for alpha and IAA concentration
y(1)=parf1TS*((SHIFTIAA+u(2))*(parff1(2)+x(2)^(2*parf1(3)))+K(1)*(parf1(1)+parf1(4)*x(2)^(2*parf1(3))))/(parf1(5)+parf1(6)*K(3)+(parf1(9)*K(3)+1)*x(2)^(2*parf1(3)))-parf1TS*x(1);
y(2)=parf2TS*((SHIFTalpha+u(1))*(parff2(2)+x(1)^(2*parf2(3)))+K(2)*(parf2(1)+parf2(4)*x(1)^(2*parf2(3))))/(parf2(5)+parf2(6)*K(4)+(parf2(9)*K(4)+1)*x(1)^(2*parf2(3)))-parf2TS*x(2);

%Defining ODEs for sensor regulation
load('ParameterEvalIAADown2xGFP')
%Adjustment from hours to minutes
parf(1:3)=parf(1:3)/60;
parf(6:end)=parf(6:end)/60;
y(3)=parf(1)*x(1)-parf(2)*x(3);
y(4)=parf(3)/(parf(4)+x(3)^parf(5))-parf(6)*x(4);
y(5)=parf(9)+parf(7)*x(4)-parf(8)*x(5);

y=y';
end