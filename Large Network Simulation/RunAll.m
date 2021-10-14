%%%Run these scripts to reproduce Figure 4 and 5

%First, run the 2-node, 3-node, and 4-node network simulations
run('Simulation/TwoNodeNetworkSimulator.m')
run('Simulation/ThreeNodeNetworkSimulator.m')
run('Simulation/FourNodeNetworkSimulator.m')

%Then, manually copy the generated files ('TwoNodeNetworkSimulation',
%'ThreeNodeNetworkSimulation', etc) to the 'LogicGatesAnalysis' folder and
%also to the 'NonMonotonicAnalysis' folder [these files are too big for me
%to plouad them!]
%Then run
run('LogicGatesAnalysis/NetworkAnalysisGates.m') %To find all the gates and optimize their values (Figure 4) 
run('NonMonotonicAnalysis/NonMonotonicFinder') %To find all the nonmonotonic behaviors and optimize them (Figure 5) 