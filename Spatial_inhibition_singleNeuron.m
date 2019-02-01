
clear;
% close all;

%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Design simulation parameters

stimRange = 1:1:30;
inhRange = [5 30 100 150]; %range of 20-40 seems like a working range?

%% Build dendritic arbor
[connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
totalCompartments = length(compartmentIDs);
shaftCmpt = compartmentIDs(1,find(compartmentIDs(2,:)~=4));
spineCmpt = compartmentIDs(1,find(compartmentIDs(2,:)==4));

%% Storage variables
noInhibition = zeros(length(spineCmpt),length(stimRange),miscParams.time);
yesInhibition = zeros(length(spineCmpt),length(shaftCmpt),length(stimRange),length(inhRange),miscParams.time);
inhDiff = zeros(length(spineCmpt),length(shaftCmpt),length(stimRange),length(inhRange));
inhDiffAOC = zeros(length(spineCmpt),length(shaftCmpt),length(stimRange),length(inhRange));
%% Run stimulation (+/- inhibition at each shaft, excitation at each spine)
tstart = tic;
for sp = 1:length(spineCmpt)
    sp
    spine = spineCmpt(sp);
    
    for stim = 1:length(stimRange)
        stimLvl = stimRange(stim);
        
        iExcite = zeros(totalCompartments,miscParams.time);
        iExcite(spine,150:200) = stimLvl;
        
        iInhibit = zeros(totalCompartments,miscParams.time);
        
        input.excitation = iExcite;
        input.inhibition = iInhibit;
        
        %run simulation
        voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
        noInhibition(sp,stim,:) = voltages(1,:);
        
        for sh = 1:length(shaftCmpt)
            shaft = shaftCmpt(sh);
            
            for inh = 1:length(inhRange)
                inhLvl = inhRange(inh);
                
                iInhibit = zeros(totalCompartments,miscParams.time);
                iInhibit(shaft,150:200) = inhLvl;
                input.inhibition = iInhibit;
                
                %run simulation
                voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
                yesInhibition(sp,sh,stim,inh,:) = voltages(1,:);
                
                %effect of inhibition is difference in max voltage of soma (positive value)
                inhDiff(sp,sh,stim,inh) = max(noInhibition(sp,stim,:)) - max(yesInhibition(sp,sh,stim,inh,:));
                inhDiffAOC(sp,sh,stim,inh) = sum(noInhibition(sp,stim,:)+60) - sum(yesInhibition(sp,sh,stim,inh,:)+60);
            end
        end
    end
    telapsed = toc(tstart)
end
g = graph(connectome);
%% Visualize data (COI = compartment of interest)
stimCOI = 40;
excLevel = 7;
inhLevel = 30;

stimIndex = find(spineCmpt==stimCOI);
excIndex = find(stimRange==excLevel);
inhIndex = find(inhRange==inhLevel);

map = jet(1000);
high = max(inhDiff(stimIndex,:,excIndex,inhIndex));
low = min(inhDiff(stimIndex,:,excIndex,inhIndex));
g = graph(connectome);
figure; h = plot(g);
for c = 1:length(shaftCmpt)
    inhCOI = shaftCmpt(c);
    difference = inhDiff(stimIndex,inhCOI,excIndex,inhIndex);
    relDiff = round(999*(difference-low)/(high-low)+1);
    cmptColor = map(relDiff,:);
    highlight(h,inhCOI);
    highlight(h,inhCOI,'NodeColor',cmptColor);
    title(['Excitation on ' num2str(stimCOI) ', intensity = ' num2str(excLevel) ', inhLevel = ' num2str(inhLevel)]);
end
highlight(h,stimCOI,'Marker','s');
highlight(h,stimCOI);
highlight(h,stimCOI);
