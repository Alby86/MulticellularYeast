%%%This code produces all the plots needed for Figure 1

clear all
close all
clc

%Simulation files
Regulation={'AlphaDownGFPsim','AlphaUpGFPsim','BetaDownGFPsim','BetaUpGFPsim','IAADown2xGFPsim','IAAUpGFPsim'};
%Data files
Datafiles={'AlphaDownGFP','AlphaUpGFP','BetaDownGFP','BetaUpGFP','IAADown2xGFP','IAAUpGFP'};

%Load background
load('BackGroundAdj','Correction');

%%%The following are 4 options on data plotting: uncomment if interested
% for j=1:4
%     figure
%     for i=1:length(Regulation)
%         load(Regulation{i})
%         load(Data{i})
%         if j==1 %Case 1: standard raw data
%         elseif j==2  %Case 2: background normalized to 1
%             Y=Y*Correction;
%             Average=Average*Correction;
%             Error=Error*Correction;
%         elseif j==3 %Case 3: background normalized to 1 and substract the minimum
%             mm=min(min(Average));
%             Y=(Y-mm)*Correction;
%             Average=(Average-mm)*Correction;
%             Error=Error*Correction;
%         elseif j==4 %Case 4: background normalized to 1 and wild type substraction
%             Y=Y*Correction-1;
%             Average=Average*Correction-1;
%             Error=Error*Correction;
%         end
%         Input(1)=0.1*Input(2);
% %         input(1)=[];
% %         Y(1)=[];
%         subplot(2,3,i)
%         errorbar(Input,Average(end,:),Error(end,:),'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
%         hold on
%         plot(input,Y,'LineWidth',2,'DisplayName','Simulation')
%         set(gca, 'XScale','log')
%         set(gca,'FontSize',15)
%         xlim([min(min(Input),min(input)) max(max(Input),max(input))])
%         ylim([0 max(max(max(Average)),max(Y))])
%         box off
%     end
% end
% 
% for i=1:length(Regulation)
%     figure
%     load(Regulation{i})
%     load(Datafiles{i})
%     Y=Y*Correction-1;
%     Average=Average*Correction-1;
%     Error=Error*Correction;
%     Input(1)=0.1*Input(2);
%     errorbar(Input,Average(end,:),Error(end,:),'MarkerSize',8,'Marker','diamond','LineWidth',3,'LineStyle','none','DisplayName','Data')
%     hold on
%     plot(input,Y,'LineWidth',3,'DisplayName','Simulation')
%     set(gca, 'XScale','log')
%     set(gca,'FontSize',15)
%     xlim([min(min(Input),min(input)) max(max(Input),max(input))])
%     ylim([0 max(max(max(Average(end,:)+Error(end,:))),max(Y))])
%     box off
%     set(gca,'XTick',[0.01 0.1 1 10 100 1000 10000]);
%     set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'})
%     savefig(Regulation{i})
%     saveas(gcf,[Regulation{i} '.svg'])
% end

%Simulation files for MISO
Regulation={'AlphaUpGFPDownIAAsim','IAAUpGFPDownAlphasim'};
%Data files
Datafiles={'AlphaUpGFPDownIAA','IAAUpGFPDownAlpha'};

for i=1:length(Regulation)
    figure
    load(Regulation{i})
    load(Datafiles{i})
    Y=Y*Correction-1;
    Data=Data*Correction-1;
    Error=Error*Correction;
    AA=unique(Input(:,2));
    for j=1:length(AA)
        c=find(AA(j)==Input(:,2));
        Input(c(Input(c,1)==0),1)=min(Input(c(Input(c,1)>0),1))*0.005;
        errorbar(Input(c,1),Data(end-1,c),Error(end-1,c),'MarkerSize',8,'Marker','diamond','LineWidth',3,'LineStyle','none','DisplayName','Data')
        hold on
        plot(input,Y(j,:),'LineWidth',3,'DisplayName','Simulation')
    end
    set(gca, 'XScale','log')
    set(gca,'FontSize',15)
    xlim([min(min(Input(:,1)),min(input)) max(max(Input(:,1)),max(input))])
    ylim([0 max(max(max(Data(end,:)+Error(end,:))),max(max(Y)))])
    box off
    set(gca,'XTick',[0.01 0.1 1 10 100 1000 10000]);
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'})
    savefig(Regulation{i})
    saveas(gcf,[Regulation{i} '.svg'])
end


%Simulation files for MISO: Time Series
Regulation={'AlphaUpGFPDownIAA','IAAUpGFPDownAlpha'};

for i=1:length(Regulation)
    figure
    load(Regulation{i})
    Data=Data*Correction-1;
    Error=Error*Correction;
    for j=1:length(Input)
        yEnd{j}=yEnd{j}*Correction-1;
        errorbar(Time(1:end-1),Data(1:end-1,j),Error(1:end-1,j),'MarkerSize',8,'Marker','diamond','LineWidth',2,'LineStyle','none','DisplayName','Data')
        hold on
        c=find(t{j}>Time(end-1),1,'first');
        plot(t{j}(1:c),yEnd{j}(1:c),'LineWidth',2,'DisplayName','Simulation')
    end
    set(gca,'FontSize',15)
    xlim([min(min(Input(:,1)),min(input)) max(max(Input(:,1)),max(input))])
    ylim([0 max(max(max(Data(end,:)+Error(end,:))),max(max(Y)))])
    box off
    set(gca,'XTick',[0.01 0.1 1 10 100 1000 10000]);
    set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'})
    savefig(Regulation{i})
    saveas(gcf,[Regulation{i} '.svg'])
end