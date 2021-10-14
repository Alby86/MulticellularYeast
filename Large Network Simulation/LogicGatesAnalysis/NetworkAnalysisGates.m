%%%This code analyze the simulated N-node networks to find logic gates

clear all
close all
clc

%This code loops over the Network datasets already simulated
SET={'Two','Three','Four'};
Nset=[2 3 4];

Tol=(0.1:0.1:1.5); %Tol+1 is the fold-change increase with respect to the minimum
MaxMax=[];
BestBest=[];

for j=1:length(SET)

    %Simulation Data Loading
    load([SET{j},'NodeNetworkSimulation'])  
    
    %Space pre-allocation
    AlphaDownGateNTAll=zeros(1,Nset(j)+3);
    AlphaUpGateNTAll=zeros(1,Nset(j)+3);
    IAADownGateNTAll=zeros(1,Nset(j)+3);
    IAAUpGateNTAll=zeros(1,Nset(j)+3);
    miso1XGateNTAll=zeros(1,Nset(j)+3);
    miso2XGateNTAll=zeros(1,Nset(j)+3);
    
    for i=1:length(Tol)        
        
        %Computing the number of gates
        [Ngates(i,:), AlphaDownGATE, IAADownGATE, AlphaUpGATE, IAAUpGATE, miso1GATE, miso2GATE]=LogicGatesFinder(GATE,GATElabel,Tol(i));
        %A typical output above is: [strain#1 strain#2...strain#Nset(j), 0/1/2 (types of input: (0=BA; 1=BI; 2=AI)), Type of gate (1=AND,2=OR,3=NAND,4=NOR,5=XOR,6=XNOR), metric value] 
        %Identifying non-trivial gates for SISO sensors
        AlphaDownGateNT=AlphaDownGATE(AlphaDownGATE(:,Nset(j)+1)==1,:);
        AlphaUpGateNT=AlphaUpGATE(AlphaUpGATE(:,Nset(j)+1)==1,:);
        IAADownGateNT=IAADownGATE(IAADownGATE(:,Nset(j)+1)==0,:);
        IAAUpGateNT=IAAUpGATE(IAAUpGATE(:,Nset(j)+1)==0,:);
        %Organizing gates across different tolerance values and eliminating
        %duplicates
        [Lia, ~]=ismember(AlphaDownGateNT, AlphaDownGateNTAll, 'rows');
        AlphaDownGateNTAll=[AlphaDownGateNTAll; AlphaDownGateNT(Lia==0,:)];
        [Lia, ~]=ismember(AlphaUpGateNT, AlphaUpGateNTAll, 'rows');
        AlphaUpGateNTAll=[AlphaUpGateNTAll; AlphaUpGateNT(Lia==0,:)];
        [Lia, ~]=ismember(IAADownGateNT, IAADownGateNTAll, 'rows');
        IAADownGateNTAll=[IAADownGateNTAll; IAADownGateNT(Lia==0,:)];
        [Lia, ~]=ismember(IAAUpGateNT, IAAUpGateNTAll, 'rows');
        IAAUpGateNTAll=[IAAUpGateNTAll; IAAUpGateNT(Lia==0,:)];
    
        %Counting the number of found gates for each value of Tol
        Ngates2(i,:)=[size(AlphaDownGateNT,1) size(AlphaUpGateNT,1) size(IAADownGateNT,1) size(IAAUpGateNT,1)];
        for k=1:6
            Ngates2type(i,k)=sum(AlphaDownGateNT(:,end)==k);
            Ngates2type(i,k+6)=sum(AlphaUpGateNT(:,end)==k);
            Ngates2type(i,k+12)=sum(IAADownGateNT(:,end)==k);
            Ngates2type(i,k+18)=sum(IAAUpGateNT(:,end)==k);
        end
    
        %For MISO systems, I consider only the XOR or XNOR gates as
        %interesting
        miso1XGateNT=miso1GATE(miso1GATE(:,end)==5 | miso1GATE(:,end)==6,:);
        miso2XGateNT=miso2GATE(miso2GATE(:,end)==5 | miso2GATE(:,end)==6,:);
        if isempty(miso1XGateNT)
        else
            [Lia, ~]=ismember(miso1XGateNT, miso1XGateNTAll, 'rows');
            miso1XGateNTAll=[miso1XGateNTAll; miso1XGateNT(Lia==0,:)];
        end
        if isempty(miso2XGateNT)
        else
            [Lia, ~]=ismember(miso2XGateNT, miso2XGateNTAll, 'rows');
            miso2XGateNTAll=[miso2XGateNTAll; miso2XGateNT(Lia==0,:)];
        end
    end
    
    %Eliminating the zeros in the first row
    AlphaDownGateNTAll(1,:)=[];
    AlphaUpGateNTAll(1,:)=[];
    IAADownGateNTAll(1,:)=[];
    IAAUpGateNTAll(1,:)=[];
    miso1XGateNTAll(1,:)=[];
    miso2XGateNTAll(1,:)=[];
    
    %Now filtering all the gates that don't show at least 1.2 fold change
    %between the minimum and the maximum according to our metric
    FoldTol=1;
    IAAUpGateNTAll(IAAUpGateNTAll(:,end)<FoldTol,:)=[];
    IAADownGateNTAll(IAADownGateNTAll(:,end)<FoldTol,:)=[];
    AlphaUpGateNTAll(AlphaUpGateNTAll(:,end)<FoldTol,:)=[];
    AlphaDownGateNTAll(AlphaDownGateNTAll(:,end)<FoldTol,:)=[];
    miso1XGateNTAll(miso1XGateNTAll(:,end)<FoldTol,:)=[];
    miso2XGateNTAll(miso2XGateNTAll(:,end)<FoldTol,:)=[];
    
    
    AlphaDownSize=size(AlphaDownGateNTAll,1);
    AlphaUpSize=size(AlphaUpGateNTAll,1);
    IAADownSize=size(IAADownGateNTAll,1);
    IAAUpSize=size(IAAUpGateNTAll,1);
    miso1Size=size(miso1XGateNTAll,1);
    miso2Size=size(miso1XGateNTAll,1);
    AllGate=[[AlphaDownGateNTAll ones(AlphaDownSize,1)]; [IAADownGateNTAll 2*ones(IAADownSize,1)]; [AlphaUpGateNTAll 3*ones(AlphaUpSize,1)]; [IAAUpGateNTAll 4*ones(IAAUpSize,1)]; [miso1XGateNTAll 5*ones(miso1Size,1)]; [miso2XGateNTAll 6*ones(miso2Size,1)]];
    
    %Finding the best gate for each gate type over all the sensor and
    %combinations
    for i=1:4
        AllG=AllGate(AllGate(:,Nset(j)+2)==i,:);
        if isempty(AllG)
            Max(i,:)=[0 0];
            Best(i,:)=zeros(1,size(AllG,2));
        else
            [Max(i,:),IndMax]=max(AllG(:,Nset(j)+3));
            Best(i,:)=AllG(IndMax,:); 
            Max(i,2)=Best(i,end);
        end
    end
    MaxMax=[MaxMax; [Max Nset(j)*ones(size(Max,1),1)]]; %Collect the Max values and add the network node size to the matrix
    BestBest=[BestBest; [Best zeros(size(Best,1),max(Nset)+4-size(Best,2))]]; %Collect the strain combinations
        
%     filename=['LogicGates' num2str(Nset(j)) 'Node'];
%     save(filename,'Ngates','Ngates2','Ngates2type','AlphaDownGateNTAll','AlphaUpGateNTAll','IAADownGateNTAll','IAAUpGateNTAll','miso1XGateNTAll','miso2XGateNTAll','Max','Best')
    clear AlphaDownGateNTAll AlphaUpGateNTAll IAADownGateNTAll IAAUpGateNTAll miso1XGateNTAll miso2XGateNTAll
    clear AlphaDownGateNT AlphaUpGateNT IAADownGateNT IAAUpGateNT miso1XGateNT miso2XGateNT
    clear AlphaDownGate AlphaUpGate IAADownGate IAAUpGate miso1XGate miso2XGate
    clear Max Best
    
end    

%Once the gates have been identified, we automatically run optimization
%on the optimal ones
%%%Defining the gates
Logic=[0 0 0 1; %AND
       0 1 1 1; %OR
       1 1 1 0; %NAND
       1 0 0 0; %NOR
       0 1 1 0; %XOR
       1 0 0 1]; %XNOR
%Optimizing for stochiometry each gate according to its network size and topology
for i=1:size(MaxMax,1)    
    i
    if MaxMax(i,1)==0
        OptimalRatioLogic{i}=[];
        OptimalRatioLogicSym{i}=[];
    else
        [OptimalRatioLogic{i},OptimalRatioLogicSym{i}]=OptimizeGate(BestBest(i,1:(4+floor((i-1)/4))),MaxMax(i,2),Logic,MaxMax(i,3));
    end
end

filename='OptimalLogicGates';
save(filename,'OptimalRatioLogic','OptimalRatioLogicSym','MaxMax','BestBest')