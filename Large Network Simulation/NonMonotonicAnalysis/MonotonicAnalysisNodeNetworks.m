function [InterestingNet1,TotInteresting,InterestingInputs,InterestingNet1F,InterestingSeries,InterestingSeriesF]=MonotonicAnalysisNodeNetworks(Y1,Yindex,Tol,N)

%The inputs are:
% Y1 is the simulation network: rows are time points and colomns are different strain combinations
% Yindex matches each simulation row from Y1 to the correspondent strain combination (rows 1:N), and input concentration (alpha, IAA, beta)
% Tol is the tolerance to consider the behavior as interesting or not (if the fold change is lower than Tol+1, then it is not considered interesting)
% N is the number of strains selected for the network generation (2-node network, 3-node network,..., N-node network)

Y1=real(Y1);

Ydiff=diff(Y1')';
pos = Ydiff>0;
changes = xor(pos(:,1:end-1),pos(:,2:end));
num = sum(changes')';

%Filter all time series that don't have much dynamic
% Tol=3; %This indicates that interesting dynamic needs to have at least a 3% value change to be considered so
newInd=DynamicTol(Y1,num,changes,Tol); %Identify the TS with non-monotonic dynamic above the threshold

%Now we filter the networks that use the same parts but are counted more
%than once because of different inducer concentrations
[InterestingNet1,Ind1]=unique(Yindex(newInd(:,1)>0,1:N),'rows');
InterestingSeries=Y1(newInd(:,1)>0,:);
InterestingSeries=InterestingSeries(Ind1,:);
for i=1:size(InterestingNet1,1)
%     cc=find(Yindex(:,1:N)==InterestingNet1(i,1:N));
    [tf,idx] = ismember(Yindex(:,1:N),InterestingNet1(i,1:N));
    cc=find(sum(idx')==(N*(N+1))/2);
    [dd,II]=max(newInd(cc,2));
    InterestingNet1(i,N+1:N+2)=[dd,cc(II)]; %Assign to the 3rd and 4th entry of the array the values corresponding to the maximum nonmonotonic gap and its Yindex entry
    Lia=ismember(Yindex(:,1:N), InterestingNet1(i,1:N), 'rows');
    aa=newInd(:,1).*Lia;
    dd=newInd(:,2).*Lia;
    InterestingInputs{i}=[newInd(find(dd>0),2) Yindex(aa>0,N+1:size(Yindex,2))];    
end

if isempty(InterestingNet1)
    InterestingInputs=[];
    InterestingSeries=[];
end
%Now we see which of these networks use the 2 strains in parallel or in
%series, and if only one input is needed
load('ParameterListF','input','output')
ParallelSeriesNet1=ParSerDyn(InterestingNet1,input,output);
InterestingNet1=[InterestingNet1 ParallelSeriesNet1'];
   
cont=0;
if isempty(InterestingNet1)
    TotInteresting=0;
    InterestingNet1F=[];
    InterestingSeriesF=[];
else
    TotInteresting=sum(InterestingNet1(:,N+3)==1);
    for i=1:size(InterestingNet1,1)
        if InterestingNet1(i,N+3)==1
            for j=1:size(InterestingInputs{i},1)
                if sum(InterestingInputs{i}(j,:)==0)==2
                    cont=cont+1;
                    InterestingNet1F(cont,1:N)=InterestingNet1(i,1:N);
                    InterestingSeriesF(cont,:)=InterestingSeries(i,:);
                    if size(InterestingInputs{i}(j,:),2)==4 %Time series case
                        InterestingNet1F(cont,N+1:N+4)=InterestingInputs{i}(j,:);
                    elseif size(InterestingInputs{i}(j,:),2)==3 %Steady State case
                        InterestingNet1F(cont,N+1:N+3)=InterestingInputs{i}(j,:);
                    end
                end
            end
        end
    end
    if cont==0
        InterestingNet1F=[];
        InterestingSeriesF=[];
    else
        [aa,i1]=unique(InterestingNet1F(:,1:N),'rows');
        if size(InterestingNet1F,2)==N+4 %Time series case
            InterestingNet1F=[aa InterestingNet1F(i1,N+1:N+4)];            
        elseif size(InterestingNet1F,2)==N+3 %Steady State case
            InterestingNet1F=[aa InterestingNet1F(i1,N+1:N+3)];
        end
        InterestingSeriesF=InterestingSeriesF(i1,:);
    end
end

end