%%%Cost function for the curve fitting program: minimises the cost over
%%%steady state and time courses

function y=CostFunAlphaDownGfpRep3(par)

%Load data
load('AlphaDownGFP')
TS=1000;
%Run simulation to compute cost function

%     InducerTimeSeries=eval(['InducerTimeSeries' num2str(i)]); %Load the first Time Series Inducer dataset
%     TimeSeries=eval(['TimeSeries' num2str(i)]); %Load the first Time Series dataset
%     TimeSeries=resampleFN(TimeSeries);%First of all we resample the data to get it better fitted for parameter estimation
    Err1=0;
    for j=1:length(Input)
        [~,yODE1]=ode15s(@Repress,[0 TS],[0 0 0],odeset('refine',10),0,par); %Simulating the system for different input
        [t,yODE]=ode15s(@Repress,[Time(1) Time(end)],yODE1(end,:),odeset('refine',10),Input(j),par);
        y=binlin(t,yODE(:,3),Time); %Resampling to the same time point as the data
        Err1(j)=sqrt(sum((log(y)-log(DataRep3(:,j))).^2)); %Cost function contribution of the first time series
    end
%     Err=sum(Err1); %Cost function contribution of the first set of time series



y=sqrt(sum(Err1)); %Cost function evaluation

end