%%%This is the standard model for gene activated by some signallig
%%%pathway

function y=ActiveIAAUpGFPDownAlpha(t,x,u,v,parf)

if u==0 %No IAA
    x(1)=0;
end
if v==0 %No alpha
    x(2)=0;
end

% parf(3)=1.8517e-05;

y(1)=parf(1)*u-parf(2)*x(1);
y(2)=parf(3)*v-parf(4)*x(2);
y(3)=parf(5)*x(1)^parf(6)/(parf(7)+x(1)^parf(6))-parf(8)*x(3);
y(4)=parf(9)*x(2)^parf(10)/(parf(11)+x(2)^parf(10))-parf(12)*x(4);
y(5)=(parf(13)+parf(14)*x(3))/(1+parf(15)*x(4))-parf(16)*x(5);

y=y';

end