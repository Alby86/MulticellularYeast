%%%This plot is to draw the phase portrait of the bistable switch system

clear all
close all
clc

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

%And the other parameters
K=[0.2 42.5 25 0.7]; %For K4=0.8, the 'high' state becomes unstable
SHIFTIAA=200;
SHIFTalpha=0;
x=[(0:10:100)];% (200:100:1000) (2000:1000:10000) (10000:10000:50000)];
[x1, x2] = meshgrid(logspace(-1,4,100),logspace(-1,4,100));
y = linspace(1,10000,50);
y2 = linspace(1,30,10);
[x1, x2] = meshgrid(y, y2);
% [x1, x2] = meshgrid((0:10:100),(100:1000:50000));
x1dot=parf1TS*(SHIFTIAA*(parff1(2)+x2.^(2*parf1(3)))+K(1)*(parf1(1)+parf1(4)*x2.^(2*parf1(3))))./(parf1(5)+parf1(6)*K(3)+(parf1(9)*K(3)+1)*x2.^(2*parf1(3)))-parf1TS*x1;
x2dot=parf2TS*(SHIFTalpha*(parff2(2)+x1.^(2*parf2(3)))+K(2)*(parf2(1)+parf2(4)*x1.^(2*parf2(3))))./(parf2(5)+parf2(6)*K(4)+(parf2(9)*K(4)+1)*x1.^(2*parf2(3)))-parf2TS*x2;
x1dot=x1dot./sqrt(x1dot.^2+x2dot.^2);
x2dot=x2dot./sqrt(x1dot.^2+x2dot.^2);
figure
quiver(x1/1000,x2,x1dot/20,x2dot,'MaxHeadSize',0.03,'AutoScaleFactor',0.5)
set(gca, 'YScale', 'log','XScale','log')
hold on
%Now plot the nullclines
x1=(0:0.1:10000);
x2=(0:0.1:30000);
x2=(0:0.1:100);
x1null=parf1TS*(SHIFTIAA*(parff1(2)+x2.^(2*parf1(3)))+K(1)*(parf1(1)+parf1(4)*x2.^(2*parf1(3))))./(parf1(5)+parf1(6)*K(3)+(parf1(9)*K(3)+1)*x2.^(2*parf1(3)))/parf1TS;
x2null=parf2TS*(SHIFTalpha*(parff2(2)+x1.^(2*parf2(3)))+K(2)*(parf2(1)+parf2(4)*x1.^(2*parf2(3))))./(parf2(5)+parf2(6)*K(4)+(parf2(9)*K(4)+1)*x1.^(2*parf2(3)))/parf2TS;
% figure
plot(x1/1e3,x2null)
plot(x1null/1e3,x2)

xlim([-1e-1 3])
ylim([0 30])


%Plot individual stable points with phase portrait
%First 2 equilibria (saddle and stable)
y2 = linspace(13,20,10);
y = linspace(8,20,10);
%More zoom
y = linspace(9,11.5,15);
y2 = linspace(18,21,10);
[x1, x2] = meshgrid(y, y2);
x1dot=parf1TS*(SHIFTIAA*(parff1(2)+x2.^(2*parf1(3)))+K(1)*(parf1(1)+parf1(4)*x2.^(2*parf1(3))))./(parf1(5)+parf1(6)*K(3)+(parf1(9)*K(3)+1)*x2.^(2*parf1(3)))-parf1TS*x1;
x2dot=parf2TS*(SHIFTalpha*(parff2(2)+x1.^(2*parf2(3)))+K(2)*(parf2(1)+parf2(4)*x1.^(2*parf2(3))))./(parf2(5)+parf2(6)*K(4)+(parf2(9)*K(4)+1)*x1.^(2*parf2(3)))-parf2TS*x2;
x1dot=x1dot./sqrt(x1dot.^2+x2dot.^2);
x2dot=x2dot./sqrt(x1dot.^2+x2dot.^2);
figure
quiver(x1,x2,x1dot,x2dot,'MaxHeadSize',0.2,'AutoScaleFactor',0.5)
hold on
xx1=(8:0.1:20);
xx2=(12:0.1:21);
x1null=parf1TS*(SHIFTIAA*(parff1(2)+xx2.^(2*parf1(3)))+K(1)*(parf1(1)+parf1(4)*xx2.^(2*parf1(3))))./(parf1(5)+parf1(6)*K(3)+(parf1(9)*K(3)+1)*xx2.^(2*parf1(3)))/parf1TS;
x2null=parf2TS*(SHIFTalpha*(parff2(2)+xx1.^(2*parf2(3)))+K(2)*(parf2(1)+parf2(4)*xx1.^(2*parf2(3))))./(parf2(5)+parf2(6)*K(4)+(parf2(9)*K(4)+1)*xx1.^(2*parf2(3)))/parf2TS;

plot(xx1,x2null)
plot(x1null,xx2)


%Second stable equilibria
y = linspace(1.4*1e3,1.8*1e3,20);
y2 = linspace(0.7,1.1,15);
[x1, x2] = meshgrid(y, y2);
% [x1, x2] = meshgrid((0:10:100),(100:1000:50000));
x1dot=parf1TS*(SHIFTIAA*(parff1(2)+x2.^(2*parf1(3)))+K(1)*(parf1(1)+parf1(4)*x2.^(2*parf1(3))))./(parf1(5)+parf1(6)*K(3)+(parf1(9)*K(3)+1)*x2.^(2*parf1(3)))-parf1TS*x1;
x2dot=parf2TS*(SHIFTalpha*(parff2(2)+x1.^(2*parf2(3)))+K(2)*(parf2(1)+parf2(4)*x1.^(2*parf2(3))))./(parf2(5)+parf2(6)*K(4)+(parf2(9)*K(4)+1)*x1.^(2*parf2(3)))-parf2TS*x2;
x1dot=x1dot./sqrt(x1dot.^2+x2dot.^2);
x2dot=x2dot./sqrt(x1dot.^2+x2dot.^2);
figure

quiver(x1/1000,x2,x1dot,x2dot,'MaxHeadSize',0.2,'AutoScaleFactor',0.5)
hold on

xx1=(1000:0.1:10000);
xx2=(0.7:0.1:15);
x1null=parf1TS*(SHIFTIAA*(parff1(2)+xx2.^(2*parf1(3)))+K(1)*(parf1(1)+parf1(4)*xx2.^(2*parf1(3))))./(parf1(5)+parf1(6)*K(3)+(parf1(9)*K(3)+1)*xx2.^(2*parf1(3)))/parf1TS;
x2null=parf2TS*(SHIFTalpha*(parff2(2)+xx1.^(2*parf2(3)))+K(2)*(parf2(1)+parf2(4)*xx1.^(2*parf2(3))))./(parf2(5)+parf2(6)*K(4)+(parf2(9)*K(4)+1)*xx1.^(2*parf2(3)))/parf2TS;

plot(xx1/1e3,x2null)
plot(x1null/1e3,xx2)