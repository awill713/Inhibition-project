
clear;
% close all;

repeats = 10;

%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Design simulation parameters

stimRange = 0.1:0.1:30;
% inhRange = [5 30 100 150]; %range of 20-40 seems like a working range?
inhRange = 20; %5 10 50 200

pathX = zeros(1,repeats); %holds the length of each repeat's (arbor's) path length, since it's always from 0 to 1
pathY = cell(1,1,repeats); %has to be a cell because there is a different length for each repeat (arbor). Made 3D for concatenation/plotting reasons
exciteSomaV = zeros(repeats,length(stimRange));
inhibitSomaV = cell(1,repeats); %cell, because inhibiting compartments along paths (from spine to soma) of different length. So array/matrix couldn't work

tstart = tic;
counter = 0;
counterTotal = length(stimRange)*repeats;
for r = 1:repeats

%% Build dendritic arbor
[connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
totalCompartments = length(compartmentIDs);
shaftCmpt = compartmentIDs(1,find(compartmentIDs(2,:)~=4));

g = graph(connectome);

%% Find clustered nodes of excitation
nodesInCluster = 1;

originalSpine = datasample(compartmentIDs(1,find(compartmentIDs(2,:)==4)),1); %which spine should we find the neighbors of, chosen randomly from entire neuron

distAndComp = [];
distAndComp(1,:) = distance(originalSpine,:); %distances of all compartments from originalSpine
distAndComp(2:3,:) = compartmentIDs(:,:); %second and third rows are compartment compartment ID and compartment type, respectively
distAndComp = sortrows(distAndComp',1)'; %sort the columns based on increasing distance from originalSpine (elements of first row)
distAndComp = distAndComp(:,find(distAndComp(3,:)==4)); %only allow spines, whose compartment type (row 3) == 4
clustered = distAndComp(2,1:nodesInCluster); %take the four closest compartments (compartment ID is row 2)
clustered(2,:) = distance(clustered(1,:)); %make second row of the distance from soma
clustered = sortrows(clustered',2)'; %sort based on distance from soma
clustered = clustered(1,:); %take the first row

%% Storage variables
noInhibition = zeros(length(stimRange),miscParams.time);
yesInhibition = zeros(length(shaftCmpt),length(stimRange),length(inhRange),miscParams.time);
inhDiff = zeros(length(shaftCmpt),length(stimRange),length(inhRange));
inhDiffAUC = zeros(length(shaftCmpt),length(stimRange),length(inhRange));

inhDiffpath = cell(length(stimRange),length(inhRange));

%% Run simulation
for exc = 1:length(stimRange)
    r
    stimLvl = stimRange(exc)
    
    iExcite = zeros(totalCompartments,miscParams.time);
    iExcite(clustered,150:200) = stimLvl;
    
    iInhibit = zeros(totalCompartments,miscParams.time);
    
    input.excitation = iExcite;
    input.inhibition = iInhibit;
    
    voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
    noInhibition(exc,:) = voltages(1,:);
    
    for inh = 1:length(inhRange)
        
        inhLvl = inhRange(inh);
        
        path = shortestpath(g,originalSpine,1);
        for shaft = 2:length(path) %only cycling through shafts in pathway instead of all shafts in arbor
            theShaft = path(shaft);
%         for shaft = 1:length(shaftCmpt) %if you wanted to cycle through all shafts in arbor
%             theShaft = shaftCmpt(shaft);
        
            iExcite = zeros(totalCompartments,miscParams.time);
            iExcite(clustered,150:200) = stimLvl;
            iInhibit = zeros(totalCompartments,miscParams.time);
            iInhibit(theShaft,150:200) = inhLvl;
            
            input.excitaiton = iExcite;
            input.inhibition = iInhibit;
            
            voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
            yesInhibition(theShaft,exc,inh,:) = voltages(1,:);
            
            inhDiff(theShaft,exc,inh) = max(noInhibition(exc,:)) - max(yesInhibition(theShaft,exc,inh,:));
            inhDiffAUC(theShaft,exc,inh) = sum(noInhibition(exc,:)+60) - sum(yesInhibition(theShaft,exc,inh,:)+60);
        end
    end
    counter = counter+1;
    tnow = toc(tstart);
    secondsLeft = (counterTotal * tnow / counter) - tnow;
    hours = floor(secondsLeft/3600);
    minutes = floor(mod(secondsLeft,3600)/60);
    seconds = floor(mod(secondsLeft, 60));
    display(['Estimated time remaining: ' num2str(hours) ' hr ' num2str(minutes) ' min ' num2str(seconds) ' sec']);
end

exciteSomaV(r,:) = sum(noInhibition+60,2);
inhibitSomaV{r} = sum(yesInhibition(path(2:end),:,:,:)+60,4);


%% Store the effect of inhibition along the path from spine to soma in pathY
inhLevel = inhRange; %5, 10, 50, 200
inhIndex = find(inhRange==inhLevel);
path = shortestpath(g,originalSpine,1);
tempY = zeros(length(stimRange),length(path)-1);
for ee = 1:length(stimRange)
    tempY(ee,:) = inhDiff(path(2:end),ee,inhIndex);
end
pathX(r) = length(path)-1;
pathY{r} = tempY;
end

%% Visualize results
figure;
average = mean(exciteSomaV,1);
plot(average);
hold on;
plot(average(1)*stimRange);
xlist = 1:40:length(stimRange);
xticks(xlist);
xticklabels(stimRange(xlist));
xlabel('Excitation intensity');
title('No inhibition');


%Interpolate each pathY (effect of inhibition along path from spine to
%soma), so that they are all the same length (interpResolution) so we can
%average them appropriately. Store in "compile"
interpResolution = 101;
xInterp = linspace(0,1,interpResolution);
compile = zeros(length(stimRange),interpResolution,repeats);
for r = 1:repeats
    x = linspace(0,1,pathX(r));
    for s = 1:length(stimRange)
        newY = interp1(x,pathY{r}(s,:),xInterp);
        compile(s,:,r) = newY;
    end
end

figure;
compressed = mean(compile,3); %average the effect of inhibition along the path across all repeats
rowMax = repmat(max(compressed,[],2),1,interpResolution);
normalizedCompressed = compressed./rowMax;
imagesc(normalizedCompressed);
c = colorbar;
xlabel('Path from original spine to soma');
xticks([1 interpResolution]);
xticklabels({'Distal shaft', 'Soma'});
ylabel('Excitation level');
ylist = 1:40:length(stimRange);
yticks(ylist);
yticklabels(stimRange(ylist));
ylabel(c,'Normalized (by row) difference in inhibition effect')
title(['Inhibition = ' num2str(inhRange) ', spines/shaft = ' num2str(dendParams.spinesPerCompartment) ', nodes/cluster = ' num2str(nodesInCluster)]);

figure;
plot(normalizedCompressed(:,1));
xlabel('Excitation intensity');
xlist = 1:40:length(stimRange);
xticks(xlist);
xticklabels(stimRange(xlist));
ylabel('Normalized inhibition effect');
title(['Inhibition on spine, inhibition = ' num2str(inhRange) ', spines/shaft = ' num2str(dendParams.spinesPerCompartment) ', nodes/cluster = ' num2str(nodesInCluster)]);

% save('spatial_inhibition_clusterRepeats_4nic4spc_anySpine100_data.mat')
