
clear;
close all;

%% Load parameters
[miscParams dendParams condParams ] = loadParameters;


%% Build dendritic arbor
[connectome compartmentIDs conductanceMat] = buildDendriticArbor(dendParams);
totalCompartments = size(compartmentIDs,2);


%% Specify excitatory and inhibitory input
iExcite = zeros(totalCompartments,miscParams.time);
iInhibit = zeros(totalCompartments,miscParams.time);
iExcite(1,150:200) = 1;

input.excitation = iExcite;
input.inhibition = iInhibit;

%% Run model
voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
figure;plot(voltages(1,:));
        