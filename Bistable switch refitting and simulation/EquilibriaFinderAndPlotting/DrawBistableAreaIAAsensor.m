%This code converts IAA distance in actual sensor values

load('SensorPars') %Notice: this is the same as finding the steady state values of the auxin sensor: I just grouped all the parameters to make it easier to read
load('DrawBistableArea')
load('BackGroundAdj','Correction')

for i=1:length(IAAsol)
    YSol(i,1)=parf(1)./(1+parf(3)*IAAsol(i,1).^parf(2))+parf(4);
    YSol(i,2)=parf(1)./(1+parf(3)*IAAsol(i,2).^parf(2))+parf(4);
    YSol(i,3)=YSol(i,2)-YSol(i,1);
end

%Assigning the nonzero values of the whole matrix
aa=reshape(ResultIAA',1,size(ResultIAA,1)*size(ResultIAA,2));
cont=0;
for i=1:length(aa)
    if aa(i)>0
        cont=cont+1;
        aa(i)=YSol(cont,3)*Correction-1;
    end
end

ResultIAAnew=reshape(aa,size(ResultIAA,2),size(ResultIAA,1))';

%And now plotting
figure
imagesc('XData',K3,'YData',K4,'CData',ResultIAAnew')

filename='DrawBistableAreaIAAsensor_Kpop_LessPerturb3';
save(filename,'K3','K4','YSol','ResultIAAnew')