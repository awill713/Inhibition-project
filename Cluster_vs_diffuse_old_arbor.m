%This script runs clustered vs diffuse simulation at a range of stimulus
%levels FOR THE OLD DENDRITIC ARBOR and plots the results as a heatmap

clear;
% close all;

%% Load parameters
[miscParams dendParams condParams ] = loadParameters;

%% Set up simulation parameters
stimLevels = 0:0.1:10;
nodesInCluster = 4;
clusts = zeros(length(stimLevels),nodesInCluster);
diffuses = zeros(length(stimLevels),nodesInCluster);

%% Run the simulation at each excitation level

for sl = 1:length(stimLevels)

    %hard code the setup of old dendritic arbor to be compatible with our
    %current system
clear connectome compartmentIDs conductanceMat distance
connectome=zeros(63,63);
connectome(1,2) = 1;
connectome(1,3) = 1;
connectome(3,[4 5]) = 1;

for x = 6:23
    connectome(x-2,x)=1;
end
for x = 24:43
    connectome(x-20,x) = 1;
end

for x = 44:63
    connectome(x-40,x) = 1;
end

connectome=connectome+triu(connectome,1)';

compartmentIDs(1,1:63) = 1:63;
compartmentIDs(2,1) = 0;
compartmentIDs(2,2) = 1;
compartmentIDs(2,3) = 2;
compartmentIDs(2,4:23) = 3;
compartmentIDs(2,24:63) = 4;

conductanceMat(1:23,1:23) = connectome(1:23,1:23) * dendParams.shaftConduct;
conductanceMat(1:23,24:63) = connectome(1:23,24:63) * dendParams.spineConduct;
conductanceMat(23:63,1:23) = connectome(23:63,1:23) * dendParams.spineConduct;

g = graph(connectome);
distance = distances(g);

totalCompartments = size(compartmentIDs,2);

%% Choose compartments for clustered and diffuse excitation conditions
%based on what was done in izkv5_clustersetup
clustered = [63 61 43 41];
diffuse = [63 60 42 41];

%% Set up clustered excitation
stimlevel = stimLevels(sl);
iExcite = zeros(4,totalCompartments,miscParams.time);
iInhibit = zeros(4,totalCompartments,miscParams.time);

for c = 1:4
    iExcite(c,clustered(1:c),150:200) = stimlevel;
end

maxV = zeros(2,size(iExcite,1)); %vector holding clustered (row 1) and diffuse (row 2) maximum soma voltages

%% Run model for clustered excitation
for exp = 1:size(iExcite,1)
%     exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(1,exp) = max(voltages(1,:));
end

%% Set up diffuse excitation

iExcite = zeros(4,totalCompartments,miscParams.time);
iInhibit = zeros(4,totalCompartments,miscParams.time);

for c = 1:4
    iExcite(c,diffuse(1:c),150:200) = stimlevel;
end


%% Run model for diffuse excitation
for exp = 1:size(iExcite,1)
%     exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(2,exp) = max(voltages(1,:));
end

clusts(sl,:) = maxV(1,:);
diffuses(sl,:) = maxV(2,:);

stimlevel
end


figure;colormap jet;
imagesc((clusts+60) ./ (diffuses+60));
title('Ratio of clustered and diffuse maximum voltage deflection');
xlabel('Number of stimulated compartments');
ylabel('Stimulation level');
colorbar
yticks(0:10:length(stimLevels));yticklabels(stimLevels(1:10:end));

figure;colormap jet;
imagesc(clusts - diffuses);
title('Difference between clustered and diffuse maximum voltage deflection');
xlabel('Number of stimulated compartments');
ylabel('Stimulation level');
colorbar
yticks(0:10:length(stimLevels));yticklabels(stimLevels(1:10:end));
