%%%This code is not a real code, but just a sequence of commands to ge thte
%%%EC50

%First, load the parameter file
%Then, run these lines:
EC(1)=parf1(2)/parf1(1)*(parf1(4))^(1/parf1(5));
EC(2)=parf2(2)/parf2(1)*(parf2(4))^(1/parf2(5));
EC(3)=parf3(2)/parf3(1)*(parf3(4))^(1/parf3(5));
mean(EC)
std(EC)

%Hill coefficient
Nset=[parf1(5),parf2(5),parf3(5)];
mean(Nset)
std(Nset)