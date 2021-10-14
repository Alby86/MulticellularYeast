function ParallelSeriesNet=ParSerDyn(InterestingNet,input,output)
%%%Function that identifies the networks that work in parallel or in series
ParallelSeriesNet=[];
for i=1:size(InterestingNet,1)
    %Pick all the networks that either match input/output or make
    %BAR1/alpha GH3/IAA
    if sum(input(InterestingNet(i,1),1:2).*output(InterestingNet(i,2),1:2))>0 | sum(input(InterestingNet(i,1),1:2).*output(InterestingNet(i,2),4:5))>0 | sum(input(InterestingNet(i,2),1:2).*output(InterestingNet(i,1),4:5))>0
        ParallelSeriesNet(i)=1; %The array contains 1 if the network is in series, a 0 if it is in parallel
    else
        ParallelSeriesNet(i)=0;
    end
end

end