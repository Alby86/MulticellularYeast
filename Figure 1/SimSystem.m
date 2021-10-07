%%%This code simulate the sensor strains and save the SS plot

clear all
close all
clc

%Parameter files
Parameters={'ParameterEvalAlphaDownGFP','ParameterEvalAlphaUpGFP','ParameterEvalBetaDownGFP','ParameterEvalBetaUpGFP','ParameterEvalIAADown2xGFP','ParameterEvalIAAUpGFP','ParameterEvalAlphaUpGFPDownIAA','ParameterEvalIAAUpGFPDownAlpha'};
%Data files
Datafiles={'AlphaDownGFP','AlphaUpGFP','BetaDownGFP','BetaUpGFP','IAADown2xGFP','IAAUpGFP','AlphaUpGFPDownIAA','IAAUpGFPDownAlpha'};

%Simulation parameters
TS=1000;

%Starting simulations
for i=1:length(Parameters)
    load(Parameters{i})
    load(Datafiles{i})
    input=(Input(1):0.1:Input(end))';
    if isempty(strfind(Parameters{i},'Up'))==0
        handle='Active';
        X0=zeros(1,3);
    elseif isempty(strfind(Parameters{i},'Down'))==0
        handle='Repress';
        X0=zeros(1,3);
    end
    if isempty(strfind(Parameters{i},'AlphaUpGFPDown'))==0
        handle='ActiveAlphaUpGFPDownIAA';
        X0=zeros(1,5);
        input=(min(Input(:,1)):0.1:max(Input(:,1)));
        [~,yODE1]=ode15s(str2func(handle),[0 TS],X0,odeset('refine',10),0,0,parf); %Simulating the system for different input
    	A=unique(Input(:,2));
        for j=1:length(A)
            for k=1:length(input)
                [t,yODE]=ode15s(str2func(handle),[Time(1) Time(end)],yODE1(end,:),odeset('refine',10),input(k),A(j),parf); %Simulating the system for different input    
                Y(j,k)=yODE(end,end);
            end
        end               
    elseif isempty(strfind(Parameters{i},'IAAUpGFPDown'))==0
        handle='ActiveIAAUpGFPDownAlpha';
        X0=zeros(1,5);
        input=(min(Input(:,1)):0.1:max(Input(:,1)));
        [~,yODE1]=ode15s(str2func(handle),[0 TS],X0,odeset('refine',10),0,0,parf); %Simulating the system for different input
    	A=unique(Input(:,2));
        for j=1:length(A)
            for k=1:length(input)
                [t,yODE]=ode15s(str2func(handle),[Time(1) Time(end)],yODE1(end,:),odeset('refine',10),input(k),A(j),parf); %Simulating the system for different input    
                Y(j,k)=yODE(end,end);
            end
        end
    else
        [~,yODE1]=ode15s(str2func(handle),[0 TS],X0,odeset('refine',10),zeros(size(input(1,:))),parf); %Simulating the system for different input
        for j=1:length(input)
            [t,yODE]=ode15s(str2func(handle),[Time(1) Time(end)],yODE1(end,:),odeset('refine',10),input(j,:),parf); %Simulating the system for different input    
            Y(j)=yODE(end,end);
        end
    end
    save([Datafiles{i}, 'sim'],'Input','input','Y')
end
