%%%This is the standard model for gene activated by some signallig
%%%pathway

function y=TwoStrainSimulatorRatio(t,x,inputs,par,input,output,K,T)

x(x<0)=0; %Just in case there are computational errors

%Assigning the parameters
alpha0=inputs(1);
IAA0=inputs(2);
beta=inputs(3);
parf1=par(1:10);
parf2=par(11:20);
input1=input(1:5);
input2=input(6:10);
output1=output(1:5);
output2=output(6:10);

alpha=(alpha0+K(1)*x(3)*output1(1)+K(2)*x(6)*output2(1))/(1+K(1)*x(3)*output1(4)+K(2)*x(6)*output2(4)); %updating alpha
IAA=(IAA0+K(1)*x(3)*output1(2)+K(2)*x(6)*output2(2))/(1+K(1)*x(3)*output1(5)+K(2)*x(6)*output2(5)); %updating IAA

%beta remains the same since it can neither be degraded nor created
u=alpha*input1(1)+IAA*input1(2)+beta*input1(3);
v=alpha*input2(1)+IAA*input2(2)+beta*input2(3);

if u==0
    x(1)=0;
end
%Strain #1
y(1)=parf1(1)*u-parf1(2)*x(1);
if parf1(end)==1
    y(2)=parf1(3)*x(1)^parf1(5)/(parf1(4)+x(1)^parf1(5))-parf1(6)*x(2);
else
    y(2)=parf1(3)/(parf1(4)+x(1)^parf1(5))-parf1(6)*x(2);
end
y(3)=parf1(9)+parf1(7)*x(2)-parf1(8)*x(3);

if v==0
    x(4)=0;
end
%Strain #2
y(4)=parf2(1)*v-parf2(2)*x(4);
if parf2(end)==1
    y(5)=parf2(3)*x(4)^parf2(5)/(parf2(4)+x(4)^parf2(5))-parf2(6)*x(5);
else
    y(5)=parf2(3)/(parf2(4)+x(4)^parf2(5))-parf2(6)*x(5);
end
y(6)=parf2(9)+parf2(7)*x(5)-parf2(8)*x(6);

%Now, load the sensors
%Alpha-factor sensor
load('ParameterEvalAlphaDownGFP')
if alpha==0 || T==1 %It's a transient simulation of individual strains or there is no input
    x(7)=0;
end
y(7)=parf(1)*alpha-parf(2)*x(7);
y(8)=parf(3)/(parf(4)+x(7)^parf(5))-parf(6)*x(8);
y(9)=parf(9)+parf(7)*x(8)-parf(8)*x(9);

%IAA sensor
load('ParameterEvalIAADown2xGFP')
if IAA==0 || T==1 %It's a transient simulation of individual strains or there is no input
    x(10)=0;
end
y(10)=parf(1)*IAA-parf(2)*x(10);
y(11)=parf(3)/(parf(4)+x(10)^parf(5))-parf(6)*x(11);
y(12)=parf(9)+parf(7)*x(11)-parf(8)*x(12);

%alpha-factor upregulated sensor
load('ParameterEvalAlphaupGFP')
if alpha==0 || T==1 %It's a transient simulation of individual strains or there is no input
    x(13)=0;
end
y(13)=parf(1)*alpha-parf(2)*x(13);
y(14)=parf(3)*x(13)^parf(5)/(parf(4)+x(13)^parf(5))-parf(6)*x(14);
y(15)=parf(9)+parf(7)*x(14)-parf(8)*x(15);

%IAA upregulated sensor
load('ParameterEvalIAAUpGFP')
if IAA==0 || T==1 %It's a transient simulation of individual strains or there is no input
    x(16)=0;
end
y(16)=parf(1)*IAA-parf(2)*x(16);
y(17)=parf(3)*x(16)^parf(5)/(parf(4)+x(16)^parf(5))-parf(6)*x(17);
y(18)=parf(9)+parf(7)*x(17)-parf(8)*x(18);

%MISO: IAA->GFP|-alpha
load('ParameterEvalIAAUpGFPDownAlpha')
if IAA==0 || T==1 %It's a transient simulation of individual strains or there is no input
    x(19)=0;
end
if alpha==0 || T==1 %It's a transient simulation of individual strains or there is no input
    x(20)=0;
end
y(19)=parf(1)*IAA-parf(2)*x(19);
y(20)=parf(3)*alpha-parf(4)*x(20);
y(21)=parf(5)*x(19)^parf(6)/(parf(7)+x(19)^parf(6))-parf(8)*x(21);
y(22)=parf(9)*x(20)^parf(10)/(parf(11)+x(20)^parf(10))-parf(12)*x(22);
y(23)=(parf(13)+parf(14)*x(21))/(1+parf(15)*x(22))-parf(16)*x(23);

%MISO: Alpha->GFP|-IAA
load('ParameterEvalAlphaUpGFPDownIAA')
if IAA==0 || T==1 %It's a transient simulation of individual strains or there is no input
    x(25)=0;
end
if alpha==0 || T==1 %It's a transient simulation of individual strains or there is no input
    x(24)=0;
end
y(24)=parf(1)*alpha-parf(2)*x(24);
y(25)=parf(3)*IAA-parf(4)*x(25);
y(26)=parf(5)*x(24)^parf(6)/(parf(7)+x(24)^parf(6))-parf(8)*x(26);
y(27)=parf(9)*x(25)^parf(10)/(parf(11)+x(25)^parf(10))-parf(12)*x(27);
y(28)=(parf(13)+parf(14)*x(26))/(1+parf(15)*x(27))-parf(16)*x(28);

y=y';

end

