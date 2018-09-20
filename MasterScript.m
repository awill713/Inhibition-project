%This entire script runs simulations of 1) increasing intensity on a single
%shaft (apical or basal) or spine (apical or basal) and 2) multiple
%simulations (simulationCount) of 1-4 (nodesInCluster) clustered or diffuse
%spines at increasing excitatory intensity (stimRange)

clear;
close all;

%% Load parameters
[miscParams dendParams condParams ] = loadParameters;

%% Build dendritic arbor
[connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
totalCompartments = size(compartmentIDs,2);

%%
% THIS SECTION COMPARES INCREASING EXCITATION INTENSITY, COMPARING
% SUB/SUPRALINEAR SUMMATION OF SPINES AND SHAFTS IN APICAL OR BASAL
% DENDRITES

%% Apical shaft - specify excitatory and inhibitory input
iExcite = zeros(10,totalCompartments,miscParams.time);
iInhibit = zeros(10,totalCompartments,miscParams.time);

stimCompartment = max(find(compartmentIDs(2,:)==3)); %most distal shaft
for intense = 1:10
    iExcite(intense,stimCompartment,150:200) = intense;
end

maxV = zeros(1,size(iExcite,1));

%% Apical shaft - run model
for exp = 1:size(iExcite,1)
    exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(exp) = max(voltages(1,:));
end
figure;plot(maxV);
hold on;
plot(-60+(maxV(1,1)+60).*(1:1:10));
legend('Observed','Predicted');
title('Apical shaft excitation is (sub)linear');

iExcite = zeros(10,totalCompartments,miscParams.time);
iInhibit = zeros(10,totalCompartments,miscParams.time);

%% Apical spine - specify excitatory and inhibitory input
stimCompartment = max(find(compartmentIDs(2,:)==4)); %most distal apical spine
for intense = 1:10
    iExcite(intense,stimCompartment,150:200) = intense;
end

maxV = zeros(1,size(iExcite,1));

%% Apical spine - run model
for exp = 1:size(iExcite,1)
    exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(exp) = max(voltages(1,:));
end
figure;plot(maxV);
hold on;
plot(-60+(maxV(1,1)+60).*(1:1:10));
legend('Observed','Predicted');
title('Apical spine excitation is supralinear');

%% Basal shaft - specify excitatory and inhibitory input
iExcite = zeros(10,totalCompartments,miscParams.time);
iInhibit = zeros(10,totalCompartments,miscParams.time);

stimCompartment = max(find(compartmentIDs(2,:)==1)); %most distal basal shaft
for intense = 1:10
    iExcite(intense,stimCompartment,150:200) = intense;
end

maxV = zeros(1,size(iExcite,1));

%% Basal shaft - run model
for exp = 1:size(iExcite,1)
    exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(exp) = max(voltages(1,:));
end
figure;plot(maxV);
hold on;
plot(-60+(maxV(1,1)+60).*(1:1:10));
legend('Observed','Predicted');
title('Basal shaft excitation is (sub)linear');

iExcite = zeros(10,totalCompartments,miscParams.time);
iInhibit = zeros(10,totalCompartments,miscParams.time);

%% Basal spine - specify excitatory and inhibitory input
iExcite = zeros(10,totalCompartments,miscParams.time);
iInhibit = zeros(10,totalCompartments,miscParams.time);

stimCompartment = max(find(connectome(max(find(compartmentIDs(2,:)==1)),:))); %most distal basal spine
for intense = 1:10
    iExcite(intense,stimCompartment,150:200) = intense;
end

maxV = zeros(1,size(iExcite,1));

%% Basal spine - run model
for exp = 1:size(iExcite,1)
    exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(exp) = max(voltages(1,:));
end
figure;plot(maxV);
hold on;
plot(-60+(maxV(1,1)+60).*(1:1:10));
legend('Observed','Predicted');
title('Basal spine excitation is supralinear');

iExcite = zeros(10,totalCompartments,miscParams.time);
iInhibit = zeros(10,totalCompartments,miscParams.time);
        
%%
% THE NEXT SECTION IS LOOKING AT THE INTEGRATION OF CLUSTERED VS DIFFUSE
% EXCITATION

%% Setup of multiple simulations

stimRange = 0:1:10; %range of intensities to use for stimulation
simulationCount = 10; %number of simulations for each intensity
nodesInCluster = 4; %up to how many compartments to stimulate in a cluster

clusterRange = zeros(length(stimRange),nodesInCluster); %matrix holding the average maximum voltage deflection for each stim level, clustered
diffuseRange = zeros(length(stimRange),nodesInCluster); %matrix holding the average maximum voltage deflection for each stim level, diffuse

%% Run simulations for each excitation level

for sl = 1:length(stimRange)
    
    stimLevel = stimRange(sl); %stim level to use for this batch of simulations
    
    clustMat = zeros(simulationCount,nodesInCluster); %maximum voltage deflections for every clustered simulation
    diffuseMat = zeros(simulationCount,nodesInCluster); %maximum voltage deflections for every diffuse simulation
    
    for count = 1:simulationCount
        
        % Build dendritic arbor
        [connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
        totalCompartments = size(compartmentIDs,2);
        
        
        %% Choose compartments for clustered and diffuse excitation conditions
        distAndComp = [];
        
        %find four clustered and diffuse spines
        %clustered:
%         originalSpine = max(find(compartmentIDs(2,:)==4));
        originalSpine = datasample(compartmentIDs(1,find(compartmentIDs(2,:)==4)),1); %which spine should we find the neighbors of, chosen randomly from entire neuron
        distAndComp(1,:) = distance(originalSpine,:); %distances of all compartments from originalSpine
        distAndComp(2:3,:) = compartmentIDs(:,:); %second and third rows are compartment compartment ID and compartment type, respectively
        distAndComp = sortrows(distAndComp',1)'; %sort the columns based on increasing distance from originalSpine (elements of first row)
        distAndComp = distAndComp(:,find(distAndComp(3,:)==4)); %only allow spines, whose compartment type (row 3) == 4
        clustered = distAndComp(2,1:nodesInCluster); %take the four closest compartments (compartment ID is row 2)
        clustered(2,:) = distance(clustered(1,:)); %make second row of the distance from soma
        clustered = sortrows(clustered',2)'; %sort based on distance from soma
        clustered = clustered(1,:); %take the first row
        
        %diffuse (randomly pick compartments whose distance from soma is equal to
        %those compartments in the clustered group):
        equidistances = distance(clustered,1); %distance from soma of compartments in clustered
        uniqueDist = unique(equidistances); %unique distances (could be duplicates)
        diffuse = [];
        %cycle through each unique distance and randomly choose the appropriate
        %number of compartments that are that distance from the soma
        for n = 1:length(unique(equidistances))
            dist = uniqueDist(n); %distance of interest
            compToChoose = find(distance(1,:)==dist); %all compartments of that distance from soma
            compToChoose = compToChoose(find(compartmentIDs(2,compToChoose)==4)); %only those who are spines
            chosen = datasample(compToChoose,sum(equidistances==dist),'Replace',false); %sample without replacement
            diffuse = [diffuse chosen]; %add the chosen compartment(s) to diffuse
        end
        diffuse(2,:) = distance(diffuse(1,:)); %make second row of distances from soma
        diffuse = sortrows(diffuse',2)'; %sort based on distance from soma
        diffuse = diffuse(1,:); %take the first row
        
        %% Set up clustered excitation
        
        iExcite = zeros(4,totalCompartments,miscParams.time);
        iInhibit = zeros(4,totalCompartments,miscParams.time);
        
        for c = 1:4
            iExcite(c,clustered(1:c),150:200) = stimLevel;
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
            iExcite(c,diffuse(1:c),150:200) = stimLevel;
        end
        
        
        %% Run model for diffuse excitation
        for exp = 1:size(iExcite,1)
            %     exp
            input.excitation = squeeze(iExcite(exp,:,:));
            input.inhibition = squeeze(iInhibit(exp,:,:));
            voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
            maxV(2,exp) = max(voltages(1,:));
        end
        
        clustMat(count,:) = maxV(1,:);
        diffuseMat(count,:) = maxV(2,:);
        
        display(['Simulation ' num2str(count) ' for stim level ' num2str(stimLevel)]);
        
    end
    
    %store the average deflections for this stim level in the appropriate matrix
    clusterRange(sl,:) = mean(clustMat);
    diffuseRange(sl,:) = mean(diffuseMat);
    
end

%% Plot results (either ratio or difference between clustered and diffuse for each stim level)
figure; colormap jet;
imagesc((clusterRange+60) ./ (diffuseRange+60));
title('Ratio of clustered and diffuse maximum voltage deflection');
xlabel('Number of stimulated compartments');
ylabel('Stimulation level');
colorbar

figure; colormap jet;
imagesc((clusterRange+60) - (diffuseRange+60));
title('Difference between clustered and diffuse maximum voltage deflection');
xlabel('Number of stimulated compartments');
ylabel('Stimulation level');
colorbar


%% Plot clustered vs diffuse data for single simulation (used prior to implementation of multi-simulation for loop
% figure;plot(maxV(1,:));
% hold on;
% plot(maxV(2,:))
% legend('Clustered','Diffuse');
% % clustered
% % diffuse
% g = graph(connectome);
% figure; h = plot(g);
% highlight(h,setdiff(clustered,diffuse),'NodeColor','g');
% highlight(h,setdiff(diffuse,clustered),'NodeColor','r');
% highlight(h,intersect(clustered,diffuse),'NodeColor','k');
% title('Green is clustered, red is diffuse, black is both');