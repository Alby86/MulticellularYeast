%This code plots the region of bistability

clear all
close all
clc

load('ParSScon_1Arm_V5')
parf1=parf;
clear parf
load('ParSScon_2Arm_V5')
parf2=parf;
clear parf
load('ParSScon_2Arm_V5WithExIAA')
parff2=parf;
clear parf PAR ERR
load('ParSScon_1Arm_V5WithExIAA')
parff1=parf;
clear parf PAR ERR

load('SolutionRange');
clear num1
SHIFT=200;
[c,I]=max(Kind(:,6)); %Solution chosen for further analysis
Ks=Kind(I,1:4);

%Computing the intervals around the Ks picked but fixing K2
K1=Ks(1); %amount of alpha-|IAA population
K2=Ks(2); %amount of IAA-|alpha population
K3=(0:0.01:30); %amount of alpha->GH3 population
K4=(0:0.01:1); %amount of IAA->BAR1 population

%Discretize the space of the input alpha
alphaD=[(0:0.01:1) (1.1:0.1:1000)];
contBist1=0;
for i=1:length(K3)
    i
    for k=1:length(K4)
        
         %And now loop over the alpha-factor discretization
         for j=1:length(alphaD)
                      
             IAA(j)= (SHIFT*(parff1(2)+alphaD(j)^(2*parf1(3)))+K1*(parf1(1)+parf1(4)*alphaD(j)^(2*parf1(3))))/(parf1(5)+parf1(6)*K3(i)+(parf1(9)*K3(i)+1)*alphaD(j)^(2*parf1(3)));
             alpha1(j)= K2*(parf2(1)+parf2(4)*IAA(j)^(2*parf2(3)))/(parf2(5)+parf2(6)*K4(k)+(parf2(9)*K4(k)+1)*IAA(j)^(2*parf2(3)));
         
         end
         Eq1=alphaD-alpha1; %Sign changes in this equation hints to the existance of a solution (the functions are continuous, so...)
         %Counting the number of sign changes
         pos1 = Eq1>0;
         changes1 = xor(pos1(1:end-1),pos1(2:end));
         num1(i,k) = sum(changes1);
         if sum(changes1)>1
             contBist1=contBist1+1;
             Y1New(contBist1,:)=Eq1;
             IndNew(contBist1,:)=[i,k];
             KindNew(contBist1,:)=[K3(i) K4(k)];
             a=find(changes1>0);
             try
                ResultAl(i,k)=alphaD(a(3))-alphaD(a(1));
                ResultIAA(i,k)=IAA(a(1))-IAA(a(3));
                Al(contBist1,1:2)=[alphaD(a(3)) alphaD(a(1))];
                IAAsol(contBist1,1:2)=[IAA(a(1)) IAA(a(3))];
             catch
                ResultAl(i,k)=0;
                ResultIAA(i,k)=0;
                Al(contBist1,1:2)=[alphaD(a(2)) alphaD(a(1))];
                IAAsol(contBist1,1:2)=[IAA(a(1)) IAA(a(2))];
             end
         else
             ResultAl(i,k)=0;
             ResultIAA(i,k)=0;
         end
         clear IAA alpha1 Eq1
    end
end
      
%And now plotting
figure
imagesc('XData',K3,'YData',K4,'CData',ResultAl')
figure
imagesc('XData',K3,'YData',K4,'CData',ResultIAA')

filename='DrawBistableArea';
save(filename,'K1','K2','K3','K4','ResultAl','ResultIAA','Al','IAAsol','KindNew','IndNew','Y1New')