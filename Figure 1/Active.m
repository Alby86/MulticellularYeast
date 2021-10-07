%%%This is the standard model for gene activated by some signallig
%%%pathway

function y=Active(t,x,u,par)

if u==0
    x(1)=0;
end
y(1)=par(1)*u-par(2)*x(1);
y(2)=par(3)*x(1)^par(5)/(par(4)+x(1)^par(5))-par(6)*x(2);
y(3)=par(9)+par(7)*x(2)-par(8)*x(3);

y=y';

end

