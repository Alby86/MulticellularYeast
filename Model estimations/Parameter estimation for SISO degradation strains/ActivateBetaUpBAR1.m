%%%This is the standard model for gene activated by some signallig
%%%pathway

function y=ActivateBetaUpBAR1(t,x,u,v,par)

v=v/(1+x(3));

load('ParameterEvalBetaUpGFP')

if u==0||x(1)==0
    x(1)=0;
    x(2)=0;
    x(3)=0;
end
%Sender strain
y(1)=parf(1)*u-parf(2)*x(1);
y(2)=parf(3)*x(1)^parf(5)/(parf(4)+x(1)^parf(5))-parf(6)*x(2);
y(3)=par(3)+par(1)*x(2)-par(2)*x(3);


if v==0
    x(4)=0;
    %x(5)=0;
%     x(6)=0;
end
load('ParameterEvalAlphaUpGFP')
%Sensor strain
y(4)=parf(1)*v-parf(2)*x(4);
y(5)=parf(3)*x(4)^parf(5)/(parf(4)+x(4)^parf(5))-parf(6)*x(5);
y(6)=parf(9)+parf(7)*x(5)-parf(8)*x(6);

y=y';

end

