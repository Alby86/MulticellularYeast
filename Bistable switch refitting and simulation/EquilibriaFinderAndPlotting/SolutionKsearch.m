%The aim of this code is to find the bistability region accross different
%population ratios given the models

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

%Discretize the parameter space of K1, amount of alpha-|IAA population
K1=(0.1:0.1:0.8);
%Discretize the parameter space of K2, amount of IAA-|alpha population
K2=5*(3.5:0.1:10);
%Discretize the parameter space of K3, amount of alpha->GH3 population
K3=(0.1:0.1:10);
%Discretize the parameter space of K4, amount of IAA->BAR1 population
K4=(0.1:0.1:0.7);

%Discretize the space of the input alpha
alphaD=[(0:0.1:1) (2:1000)];

%Using a numerical approach
%Loop over the different ratio Ks to find two values of alpha for which the
%function is zero
contBist1=0;%Counts the number of bistable solutions
SHIFT=200; %Exogenous IAA concentration -> Bad for the distance between equilibria
for w=1:length(K1)
    w
for i=1:length(K2)
    
    for k=1:length(K3)
        
        for q=1:length(K4)
            
            %And now loop over the alpha-factor discretization
            for j=1:length(alphaD)
                
                IAA(j)= (SHIFT*(parff1(2)+alphaD(j)^(2*parf1(3)))+K1(w)*(parf1(1)+parf1(4)*alphaD(j)^(2*parf1(3))))/(parf1(5)+parf1(6)*K3(k)+(parf1(9)*K3(k)+1)*alphaD(j)^(2*parf1(3)));
                alpha1(j)= K2(i)*(parf2(1)+parf2(4)*IAA(j)^(2*parf2(3)))/(parf2(5)+parf2(6)*K4(q)+(parf2(9)*K4(q)+1)*IAA(j)^(2*parf2(3)));
                
            end
            Eq1=alphaD-alpha1; %Sign changes in this equation hints to the existance of a solution (the functions are continuous, so...)
            %Counting the number of sign changes
            pos1 = Eq1>0;
            changes1 = xor(pos1(1:end-1),pos1(2:end));
            num1(w,i,k,q) = sum(changes1);
            if sum(changes1)>1
                contBist1=contBist1+1;
                Y1(contBist1,:)=Eq1;
                Ind1(contBist1,:)=[w,i,k,q];
                Kind(contBist1,1:4)=[K1(w) K2(i) K3(k) K4(q)];
                a=find(changes1>0);
             try
                Kind(contBist1,5)=alphaD(a(3))-alphaD(a(1));
                Kind(contBist1,6)=IAA(a(1))-IAA(a(3));
             catch
                Kind(contBist1,5)=alphaD(a(2))-alphaD(a(1));
                Kind(contBist1,6)=IAA(a(1))-IAA(a(2));
             end
            end
        end
    end
end
end

filename='SolutionRange';
save(filename,'num1','contBist1','Y1','Ind1','Kind')%,'MaxParRatio')%,'Best')