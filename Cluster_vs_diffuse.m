%This script runs multiple simulations (simulationCount) of 1-4
%(nodesInCluster) clustered or diffuse spines at increasing excitatory
%intensity (stimRange), and compares the clustered vs diffuse data as a
%heat map



clear;
% close all;

%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Setup of multiple simulations

stimRange = 1:0.1:10; %range of intensities to use for stimulation
simulationCount = 10; %number of simulations for each intensity
nodesInCluster = 4; %up to how many compartments to stimulate in a cluster

stimOn = 150;
stimOff = 200;

clusterRange = zeros(length(stimRange),nodesInCluster); %matrix holding the average maximum voltage deflection for each stim level, clustered
diffuseRange = zeros(length(stimRange),nodesInCluster); %matrix holding the average maximum voltage deflection for each stim level, diffuse

clusterRangeAUC = zeros(length(stimRange),nodesInCluster);
diffuseRangeAUC = zeros(length(stimRange),nodesInCluster);

%% Run simulations for each excitation level

for sl = 1:length(stimRange)
    
    stimLevel = stimRange(sl) %stim level to use for this batch of simulations
    
    clustMat = zeros(simulationCount,nodesInCluster); %maximum voltage deflections for every clustered simulation
    diffuseMat = zeros(simulationCount,nodesInCluster); %maximum voltage deflections for every diffuse simulation
    
    clustMatAUC = zeros(simulationCount,nodesInCluster);
    diffuseMatAUC = zeros(simulationCount,nodesInCluster);
    
    for count = 1:simulationCount
        
        % Build dendritic arbor
        [connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
        totalCompartments = size(compartmentIDs,2);
        
        
        %% Choose compartments for clustered and diffuse excitation conditions
        distAndComp = [];
        
        %find the clustered and diffuse spines
        %clustered:
%         originalSpine = max(find(compartmentIDs(2,:)==4)); %most distal spine
        originalSpine = datasample(compartmentIDs(1,find(compartmentIDs(2,:)==4)),1); %randomly choose spine
        
        distAndComp(1,:) = distance(originalSpine,:); %distances of all compartments from originalSpine
        distAndComp(2:3,:) = compartmentIDs(:,:); %second and third rows are compartment compartment ID and compartment type, respectively
        distAndComp = sortrows(distAndComp',1)'; %sort the columns based on increasing distance from originalSpine (elements of first row)
        distAndComp = distAndComp(:,find(distAndComp(3,:)==4)); %only allow spines, compartment type (row 3) == 4
        clustered = distAndComp(2,1:nodesInCluster); %take the "nodesInCluster" number of closest compartments (compartment ID is row 2)
        clustered(2,:) = distance(clustered(1,:)); %make second row the distance from soma
        clustered = sortrows(clustered',2)'; %sort based on distance from soma
        clustered = clustered(1,:); %take the first row, aka the spine compartments themselves
        
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
        
        iExcite = zeros(nodesInCluster,totalCompartments,miscParams.time);
        iInhibit = zeros(nodesInCluster,totalCompartments,miscParams.time);
        
        for c = 1:nodesInCluster
            iExcite(c,clustered(1:c),stimOn:stimOff) = stimLevel;
        end
        
        maxV = zeros(2,size(iExcite,1)); %vector holding clustered (row 1) and diffuse (row 2) maximum soma voltages
        somaAUC = zeros(2,size(iExcite,1));
        
        %% Run model for clustered excitation
        for exp = 1:size(iExcite,1)
            input.excitation = squeeze(iExcite(exp,:,:));
            input.inhibition = squeeze(iInhibit(exp,:,:));
            voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
            maxV(1,exp) = max(voltages(1,:));
            somaAUC(1,exp) = sum(voltages(1,150:250)-condParams.vRest);
        end
        
        %% Set up diffuse excitation
        
        iExcite = zeros(nodesInCluster,totalCompartments,miscParams.time);
        iInhibit = zeros(nodesInCluster,totalCompartments,miscParams.time);
        
        for c = 1:nodesInCluster
            iExcite(c,diffuse(1:c),stimOn:stimOff) = stimLevel;
        end
        
        
        %% Run model for diffuse excitation
        for exp = 1:size(iExcite,1)
            input.excitation = squeeze(iExcite(exp,:,:));
            input.inhibition = squeeze(iInhibit(exp,:,:));
            voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
            maxV(2,exp) = max(voltages(1,:));
            somaAUC(2,exp) = sum(voltages(1,150:250)-condParams.vRest);
        end
        
        clustMat(count,:) = maxV(1,:);
        diffuseMat(count,:) = maxV(2,:);
        
        clustMatAUC(count,:) = somaAUC(1,:);
        diffuseMatAUC(count,:) = somaAUC(2,:);
        
        display(['Simulation ' num2str(count) ' for stim level ' num2str(stimLevel)]);
        
    end
    
    %store the average deflections for this stim level in the appropriate matrix
    clusterRange(sl,:) = mean(clustMat);
    diffuseRange(sl,:) = mean(diffuseMat);
    
    clusterRangeAUC(sl,:) = mean(clustMatAUC);
    diffuseRangeAUC(sl,:) = mean(diffuseMatAUC);

    
end

%% Plot results (either ratio or difference between clustered and diffuse for each stim level)
f1 = figure; colormap jet;
imagesc(1:nodesInCluster,stimRange,(clusterRange+60) ./ (diffuseRange+60));
title('Ratio of clustered and diffuse maximum voltage deflection');
xlabel('Number of stimulated compartments');
ylabel('Stimulation level');
c = colorbar;
ylabel(c,'Ratio of clustered to diffuse maximum voltage deflection')

f2 = figure; colormap jet;
imagesc(1:nodesInCluster,stimRange,(clusterRange+60) - (diffuseRange+60));
title('Difference between clustered and diffuse maximum voltage deflection');
xlabel('Number of stimulated compartments');
ylabel('Stimulation level');
c = colorbar;
ylabel(c,'Difference between clustered and diffuse maximum voltage deflection')

f3 = figure; colormap jet;
imagesc(1:nodesInCluster,stimRange,clusterRangeAUC ./ diffuseRangeAUC);
title('Ratio of clustered and diffuse maximum AUC');
xlabel('Number of stimulated compartments');
ylabel('Stimulation level');
c = colorbar;
ylabel(c,'Ratio of clustered to diffuse maximum AUC')

f4 = figure; colormap jet;
imagesc(1:nodesInCluster,stimRange,clusterRangeAUC - diffuseRangeAUC);
title('Difference between clustered and diffuse AUC');
xlabel('Number of stimulated compartments');
ylabel('Stimulation level');
c = colorbar;
ylabel(c,'Difference between clustered and diffuse AUC')

predictedAUCclust = zeros(1,nodesInCluster);
predictedAUCdiff = zeros(1,nodesInCluster);
AUCclustMaxV = zeros(1,nodesInCluster);
AUCdiffMaxV = zeros(1,nodesInCluster);
AUCclustAUC = zeros(1,nodesInCluster);
AUCdiffAUC = zeros(1,nodesInCluster);
for i = 1:nodesInCluster
    predictedAUCclust(i) = sum((-60+(clusterRange(1,i)+60)*stimRange) - condParams.vRest);
    predictedAUCdiff(i) = sum((-60+(diffuseRange(1,i)+60)*stimRange) - condParams.vRest);
    AUCclustMaxV(i) = (sum(clusterRange(:,i) - condParams.vRest) - predictedAUCclust(i)) / predictedAUCclust(i);
    AUCdiffMaxV(i) = (sum(diffuseRange(:,i) - condParams.vRest) - predictedAUCdiff(i)) / predictedAUCdiff(i);
    
    predictedAUCclust(i) = sum(clusterRangeAUC(1,i)*stimRange);
    predictedAUCdiff(i) = sum(diffuseRangeAUC(1,i)*stimRange);
    AUCclustAUC(i) = (sum(clusterRangeAUC(:,i)) - predictedAUCclust(i)) / predictedAUCclust(i);
    AUCdiffAUC(i) = (sum(diffuseRangeAUC(:,i)) - predictedAUCdiff(i)) / predictedAUCdiff(i);
        
    if i==nodesInCluster
        f5 = figure;hold on;
        plot(stimRange,clusterRange(:,i));
        plot(stimRange,-60+(clusterRange(1,i)+60)*stimRange);
        plot(stimRange,diffuseRange(:,i));
        plot(stimRange,-60+(diffuseRange(1,i)+60)*stimRange);
        legend({'Clust obs','Clust pred','Diff obs','Diff pred'});
        title(['Somatic maximum voltage, stimulating ' num2str(nodesInCluster) ' spines']);
        xlabel('Stimulation level');
        ylabel('Maximum somatic voltage (mV)');
        
        f6 = figure;hold on;
        plot(stimRange,clusterRangeAUC(:,i));
        plot(stimRange,clusterRangeAUC(1,i)*stimRange);
        plot(stimRange,diffuseRangeAUC(:,i));
        plot(stimRange,diffuseRangeAUC(1,i)*stimRange);
        legend({'Clust obs','Clust pred','Diff obs','Diff pred'});
        title(['Somatic voltage area under curve, stimulating ' num2str(nodesInCluster) ' spines']);
        xlabel('Stimulation level');
        ylabel('mVolt * time ...?');
    end
end

f7 = figure;
subplot(2,1,1);
bar([AUCdiffMaxV' AUCclustMaxV']);
title('Nonlinearity indices of clustered and diffuse, max voltage');
ylabel('Nonlinearity index');
legend({'Diffuse','Clustered'});
subplot(2,1,2);
bar([(AUCclustMaxV - AUCdiffMaxV)']);
ylabel('Difference in nonlinearity index');
xlabel('Number of spines stimulated');

f8 = figure;
subplot(2,1,1);
bar([AUCdiffAUC' AUCclustAUC']);
title('Nonlinearity indices of clustered and diffuse, area under curve');
legend({'Diffuse','Clustered'});
ylabel('Nonlinearity index');
subplot(2,1,2);
bar([(AUCclustAUC - AUCdiffAUC)']);
ylabel('Difference in nonlinearity index');
xlabel('Number of spines stimulated');


%% Save everything

data.dendParams = dendParams;
data.miscParams = miscParams;
data.condParams = condParams;
data.stimRange = stimRange;
data.simulationCount = simulationCount;
data.nodesInCluster = nodesInCluster;
data.clusterData = clusterRange;
data.diffuseData = diffuseRange;
data.clusterDataAUC = clusterRangeAUC;
data.diffuseDataAUC = diffuseRangeAUC;
% 
% currentFolder = pwd;
% save([pwd '/Figures/cluster_vs_diffuse_data.mat'],'data');
% saveas(f1,[pwd '/Figures/clustered vs diffuse ratio.fig']);
% saveas(f2,[pwd '/Figures/clustered vs diffuse difference.fig']);
% saveas(f3,[pwd '/Figures/clustered vs diffuse AUC ratio.fig']);
% saveas(f4,[pwd '/Figures/clustered vs diffuse AUC difference.fig']);
% saveas(f5,[pwd '/Figures/clustered vs diffuse vs predicted.fig']);
% saveas(f6,[pwd '/Figures/clustered vs diffuse vs predicted AUC.fig']);
% saveas(f7,[pwd '/Figures/clustered vs diffuse nonlinearity index.fig']);
% saveas(f8,[pwd '/Figures/clustered vs diffuse AUC nonlinearity index.fig']);