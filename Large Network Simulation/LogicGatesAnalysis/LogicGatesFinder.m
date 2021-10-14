%%%This code analyse the data simulated in the N-node network and
%%%identifies strain compositions that bring about 'logic'

function [Ngates, AlphaDownGATE, IAADownGATE, AlphaUpGATE, IAAUpGATE, miso1GATE, miso2GATE]=LogicGatesFinder(GATE,GATElabel,Tol)

%%%First, let's reshape the matrices
GATE=real(GATE);

GATEAlphaDown=GATE(:,1:6:end);
GATEAlphaDown(:,8)=[];
GATEAlphaDownBA=GATEAlphaDown;
GATEAlphaDownBA(:,5:7)=[];
GATEAlphaDownBI=GATEAlphaDown;
GATEAlphaDownBI(:,7)=[];
GATEAlphaDownBI(:,3:4)=[];
GATEAlphaDownAI=GATEAlphaDown;
GATEAlphaDownAI(:,6)=[];
GATEAlphaDownAI(:,4)=[];
GATEAlphaDownAI(:,2)=[];

GATEAlphaDownAll=[GATEAlphaDownBA; GATEAlphaDownBI; GATEAlphaDownAI];

GATEIAADown=GATE(:,2:6:end);
GATEIAADown(:,8)=[];
GATEIAADownBA=GATEIAADown;
GATEIAADownBA(:,5:7)=[];
GATEIAADownBI=GATEIAADown;
GATEIAADownBI(:,7)=[];
GATEIAADownBI(:,3:4)=[];
GATEIAADownAI=GATEIAADown;
GATEIAADownAI(:,6)=[];
GATEIAADownAI(:,4)=[];
GATEIAADownAI(:,2)=[];

GATEIAADownAll=[GATEIAADownBA; GATEIAADownBI; GATEIAADownAI];

GATEAlphaUp=GATE(:,3:6:end);
GATEAlphaUp(:,8)=[];
GATEAlphaUpBA=GATEAlphaUp;
GATEAlphaUpBA(:,5:7)=[];
GATEAlphaUpBI=GATEAlphaUp;
GATEAlphaUpBI(:,7)=[];
GATEAlphaUpBI(:,3:4)=[];
GATEAlphaUpAI=GATEAlphaUp;
GATEAlphaUpAI(:,6)=[];
GATEAlphaUpAI(:,4)=[];
GATEAlphaUpAI(:,2)=[];

GATEAlphaUpAll=[GATEAlphaUpBA; GATEAlphaUpBI; GATEAlphaUpAI];

GATEIAAUp=GATE(:,4:6:end);
GATEIAAUp(:,8)=[];
GATEIAAUpBA=GATEIAAUp;
GATEIAAUpBA(:,5:7)=[];
GATEIAAUpBI=GATEIAAUp;
GATEIAAUpBI(:,7)=[];
GATEIAAUpBI(:,3:4)=[];
GATEIAAUpAI=GATEIAAUp;
GATEIAAUpAI(:,6)=[];
GATEIAAUpAI(:,4)=[];
GATEIAAUpAI(:,2)=[];

GATEIAAUpAll=[GATEIAAUpBA; GATEIAAUpBI; GATEIAAUpAI];

GATEmiso1=GATE(:,5:6:end);
GATEmiso1(:,8)=[];
GATEmiso1BA=GATEmiso1;
GATEmiso1BA(:,5:7)=[];
GATEmiso1BI=GATEmiso1;
GATEmiso1BI(:,7)=[];
GATEmiso1BI(:,3:4)=[];
GATEmiso1AI=GATEmiso1;
GATEmiso1AI(:,6)=[];
GATEmiso1AI(:,4)=[];
GATEmiso1AI(:,2)=[];

GATEmiso1All=[GATEmiso1BA; GATEmiso1BI; GATEmiso1AI];

GATEmiso2=GATE(:,6:6:end);
GATEmiso2(:,8)=[];
GATEmiso2BA=GATEmiso2;
GATEmiso2BA(:,5:7)=[];
GATEmiso2BI=GATEmiso2;
GATEmiso2BI(:,7)=[];
GATEmiso2BI(:,3:4)=[];
GATEmiso2AI=GATEmiso2;
GATEmiso2AI(:,6)=[];
GATEmiso2AI(:,4)=[];
GATEmiso2AI(:,2)=[];

GATEmiso2All=[GATEmiso2BA; GATEmiso2BI; GATEmiso2AI];

%Saving the original number to compute the metric later on
GATEAlphaDownOld=GATEAlphaDown;
GATEAlphaDownAllOld=GATEAlphaDownAll;
GATEIAADownOld=GATEIAADown;
GATEIAADownAllOld=GATEIAADownAll;
GATEAlphaUpOld=GATEAlphaUp;
GATEAlphaUpAllOld=GATEAlphaUpAll;
GATEIAAUpOld=GATEIAAUp;
GATEIAAUpAllOld=GATEIAAUpAll;
GATEmiso1Old=GATEmiso1;
GATEmiso1AllOld=GATEmiso1All;
GATEmiso2Old=GATEmiso2;
GATEmiso2AllOld=GATEmiso2All;

%Next, I normalize the values so that the minimum is 0 and the other
%entries are measured as fold changes of the minimum
for i=1:length(GATE)
    
    GATEAlphaDown(i,:)=GATEAlphaDown(i,:)/min(GATEAlphaDown(i,:))-1;
    GATEIAADown(i,:)=GATEIAADown(i,:)/min(GATEIAADown(i,:))-1;
    GATEAlphaUp(i,:)=GATEAlphaUp(i,:)/min(GATEAlphaUp(i,:))-1;
    GATEIAAUp(i,:)=GATEIAAUp(i,:)/min(GATEIAAUp(i,:))-1;
    GATEmiso1(i,:)=GATEmiso1(i,:)/min(GATEmiso1(i,:))-1;
    GATEmiso2(i,:)=GATEmiso2(i,:)/min(GATEmiso2(i,:))-1;
    
    for j=1:3
        
        GATEAlphaUpAll(i+(j-1)*length(GATE),:)=GATEAlphaUpAll(i+(j-1)*length(GATE),:)/min(GATEAlphaUpAll(i+(j-1)*length(GATE),:))-1;
        GATEIAAUpAll(i+(j-1)*length(GATE),:)=GATEIAAUpAll(i+(j-1)*length(GATE),:)/min(GATEIAAUpAll(i+(j-1)*length(GATE),:))-1;
        GATEAlphaDownAll(i+(j-1)*length(GATE),:)=GATEAlphaDownAll(i+(j-1)*length(GATE),:)/min(GATEAlphaDownAll(i+(j-1)*length(GATE),:))-1;
        GATEIAADownAll(i+(j-1)*length(GATE),:)=GATEIAADownAll(i+(j-1)*length(GATE),:)/min(GATEIAADownAll(i+(j-1)*length(GATE),:))-1;
        GATEmiso1All(i+(j-1)*length(GATE),:)=GATEmiso1All(i+(j-1)*length(GATE),:)/min(GATEmiso1All(i+(j-1)*length(GATE),:))-1;
        GATEmiso2All(i+(j-1)*length(GATE),:)=GATEmiso2All(i+(j-1)*length(GATE),:)/min(GATEmiso2All(i+(j-1)*length(GATE),:))-1;
        
    end
    
end


%Final 'logic' normalization. NOTE: this is a naive approach. Another
%possible approach is to consider if the mean of the off states is below
%the threshold but that's hard to generalize
GATEAlphaDown=(GATEAlphaDown>Tol);
GATEAlphaDownAll=(GATEAlphaDownAll>Tol);

GATEIAADown=(GATEIAADown>Tol);
GATEIAADownAll=(GATEIAADownAll>Tol);

GATEAlphaUp=(GATEAlphaUp>Tol);
GATEAlphaUpAll=(GATEAlphaUpAll>Tol);


GATEIAAUp=(GATEIAAUp>Tol);
GATEIAAUpAll=(GATEIAAUpAll>Tol);

GATEmiso1=(GATEmiso1>Tol);
GATEmiso1All=(GATEmiso1All>Tol);

GATEmiso2=(GATEmiso2>Tol);
GATEmiso2All=(GATEmiso2All>Tol);

%%%Defining the gates
Logic=[0 0 0 1; %AND
       0 1 1 1; %OR
       1 1 1 0; %NAND
       1 0 0 0; %NOR
       0 1 1 0; %XOR
       1 0 0 1]; %XNOR

%And now, comparing to find the gates

[Lia, Locb]=ismember(GATEAlphaDownAll, Logic, 'rows'); %Identifying gates using alpha-|GFP sensor
 %Computing the metric Max/min
Y=GATEAlphaDownAllOld(Lia>0,:);
Logic0=Logic(Locb(Locb>0,:),:);
M=Y.*Logic0;
m=Y.*(Logic0==0);
for i=1:size(M,1) %The metric is computed as the division of the minimum between the '1' values and the maximum between the '0' values
    a(i)=min(M(i,M(i,:)>0));
    b(i)=max(m(i,:));
end
%Add the gate information and strain info in one matrix
ind=mod(find(Lia>0),length(GATE))+length(GATE)*double(mod(find(Lia>0),length(GATE))==0);
AlphaDownGATE=[GATElabel(ind,:) floor(find(Lia>0)/length(GATE)) Locb(Lia>0)]; %Spilling out the gates in terms of [strain combo; inputs (0=BA; 1=BI; 2=AI); type of gate]
if isempty(M)
    S=[];
    AlphaDownGATE=double.empty(0,size(AlphaDownGATE,2)+1);
else
    S=a./b;
    AlphaDownGATE=[AlphaDownGATE S']; %Adding the values of the metric
end
NAlphaDownGATE=length(find(Lia>0));
clear a b S

[Lia, Locb]=ismember(GATEAlphaUpAll, Logic, 'rows'); %Identifying gates using alpha-|GFP sensor
 %Computing the metric Max/min
Y=GATEAlphaUpAllOld(Lia>0,:);
Logic0=Logic(Locb(Locb>0,:),:);
M=Y.*Logic0;
m=Y.*(Logic0==0);
for i=1:size(M,1) %The metric is computed as the division of the minimum between the '1' values and the maximum between the '0' values
    a(i)=min(M(i,M(i,:)>0));
    b(i)=max(m(i,:));
end
%Add the gate information and strain info in one matrix
ind=mod(find(Lia>0),length(GATE))+length(GATE)*double(mod(find(Lia>0),length(GATE))==0);
AlphaUpGATE=[GATElabel(ind,:) floor(find(Lia>0)/length(GATE)) Locb(Lia>0)]; %Spilling out the gates in terms of [strain combo; inputs (0=BA; 1=BI; 2=AI); type of gate]
if isempty(M)
    S=[];
    AlphaUpGATE=double.empty(0,size(AlphaUpGATE,2)+1);
else
    S=a./b;
    AlphaUpGATE=[AlphaUpGATE S']; %Adding the values of the metric
end
NAlphaUpGATE=length(find(Lia>0));
clear a b S

[Lia, Locb]=ismember(GATEIAADownAll, Logic, 'rows'); %Identifying gates using alpha-|GFP sensor
 %Computing the metric Max/min
Y=GATEIAADownAllOld(Lia>0,:);
Logic0=Logic(Locb(Locb>0,:),:);
M=Y.*Logic0;
m=Y.*(Logic0==0);
for i=1:size(M,1) %The metric is computed as the division of the minimum between the '1' values and the maximum between the '0' values
    a(i)=min(M(i,M(i,:)>0));
    b(i)=max(m(i,:));
end
%Add the gate information and strain info in one matrix
ind=mod(find(Lia>0),length(GATE))+length(GATE)*double(mod(find(Lia>0),length(GATE))==0);    
IAADownGATE=[GATElabel(ind,:) floor(find(Lia>0)/length(GATE)) Locb(Lia>0)]; %Spilling out the gates in terms of [strain combo; inputs (0=BA; 1=BI; 2=AI); type of gate]
if isempty(M)
    S=[];
    IAADownGATE=double.empty(0,size(IAADownGATE,2)+1);
else
    S=a./b;
    IAADownGATE=[IAADownGATE S']; %Adding the values of the metric
end
NIAADownGATE=length(find(Lia>0));
clear a b S

[Lia, Locb]=ismember(GATEIAAUpAll, Logic, 'rows'); %Identifying gates using alpha-|GFP sensor
 %Computing the metric Max/min
Y=GATEIAAUpAllOld(Lia>0,:);
Logic0=Logic(Locb(Locb>0,:),:);
M=Y.*Logic0;
m=Y.*(Logic0==0);
for i=1:size(M,1) %The metric is computed as the division of the minimum between the '1' values and the maximum between the '0' values
    a(i)=min(M(i,M(i,:)>0));
    b(i)=max(m(i,:));
end
%Add the gate information and strain info in one matrix
ind=mod(find(Lia>0),length(GATE))+length(GATE)*double(mod(find(Lia>0),length(GATE))==0);
IAAUpGATE=[GATElabel(ind,:) floor(find(Lia>0)/length(GATE)) Locb(Lia>0)]; %Spilling out the gates in terms of [strain combo; inputs (0=BA; 1=BI; 2=AI); type of gate]
if isempty(M)
    S=[];
    IAAUpGATE=double.empty(0,size(IAAUpGATE,2)+1);
else
    S=a./b;
    IAAUpGATE=[IAAUpGATE S']; %Adding the values of the metric
end
NIAAUpGATE=length(find(Lia>0));
clear a b S

[Lia, Locb]=ismember(GATEmiso1All, Logic, 'rows'); %Identifying gates using alpha-|GFP sensor
 %Computing the metric Max/min
Y=GATEmiso1AllOld(Lia>0,:);
Logic0=Logic(Locb(Locb>0,:),:);
M=Y.*Logic0;
m=Y.*(Logic0==0);
for i=1:size(M,1) %The metric is computed as the division of the minimum between the '1' values and the maximum between the '0' values
    a(i)=min(M(i,M(i,:)>0));
    b(i)=max(m(i,:));
end
%Add the gate information and strain info in one matrix
ind=mod(find(Lia>0),length(GATE))+length(GATE)*double(mod(find(Lia>0),length(GATE))==0);
miso1GATE=[GATElabel(ind,:) floor(find(Lia>0)/length(GATE)) Locb(Lia>0)]; %Spilling out the gates in terms of [strain combo; inputs (0=BA; 1=BI; 2=AI); type of gate]
if isempty(M)
    S=[];
    miso1GATE=double.empty(0,size(miso1GATE,2)+1);
else
    S=a./b;
    miso1GATE=[miso1GATE S']; %Adding the values of the metric
end
Nmiso1GATE=length(find(Lia>0));
clear a b S

[Lia, Locb]=ismember(GATEmiso2All, Logic, 'rows'); %Identifying gates using alpha-|GFP sensor
 %Computing the metric Max/min
Y=GATEmiso2AllOld(Lia>0,:);
Logic0=Logic(Locb(Locb>0,:),:);
M=Y.*Logic0;
m=Y.*(Logic0==0);
for i=1:size(M,1) %The metric is computed as the division of the minimum between the '1' values and the maximum between the '0' values
    a(i)=min(M(i,M(i,:)>0));
    b(i)=max(m(i,:));
end
%Add the gate information and strain info in one matrix
ind=mod(find(Lia>0),length(GATE))+length(GATE)*double(mod(find(Lia>0),length(GATE))==0);
miso2GATE=[GATElabel(ind,:) floor(find(Lia>0)/length(GATE)) Locb(Lia>0)]; %Spilling out the gates in terms of [strain combo; inputs (0=BA; 1=BI; 2=AI); type of gate]
if isempty(M)
    S=[];
    miso2GATE=double.empty(0,size(miso2GATE,2)+1);
else
    S=a./b;
    miso2GATE=[miso2GATE S']; %Adding the values of the metric
end
Nmiso2GATE=length(find(Lia>0));
clear a b S

Ngates=[NAlphaDownGATE NIAADownGATE NAlphaUpGATE NIAAUpGATE Nmiso1GATE Nmiso2GATE];

end