
%This script looks at the magnitude/effect of inhibition with various
%delays relative to the onset of excitation.


clear;
% close all;

repeats = 10;

%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Design simulation parameters

stimRange = 0.1:0.1:10;
inhRange = 20; %5 10 50 200
inhOffset = -100:10:100; %in tenths of ms (so actually -10:10 ms)

inhEffect = zeros(repeats,length(stimRange),length(inhOffset));
inhEffectAUC = zeros(repeats,length(stimRange),length(inhOffset));

%% Run simulation
tstart = tic;
for r = 1:repeats
    
%% Build dendritic arbor
[connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
totalCompartments = length(compartmentIDs);
shaftCmpt = compartmentIDs(1,find(compartmentIDs(2,:)~=4));
spineCmpt = compartmentIDs(1,find(compartmentIDs(2,:)==4));

spine = datasample(spineCmpt,1);
shaft = find(connectome(spine,:)==1);


%% Storage variables
noInhibition = zeros(length(stimRange),miscParams.time);
yesInhibition = zeros(length(stimRange),length(inhOffset),miscParams.time);
inhDiff = zeros(length(stimRange),length(inhOffset));
inhDiffAUC = zeros(length(stimRange),length(inhOffset));

%% Run simulation
for exc = 1:length(stimRange)
    stimLvl = stimRange(exc);
    
    iExcite = zeros(totalCompartments,miscParams.time);
    iExcite(spine,250:300) = stimLvl;
    
    iInhibit = zeros(totalCompartments,miscParams.time);
    
    input.excitation = iExcite;
    input.inhibition = iInhibit;
    
    voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
    noInhibition(exc,:) = voltages(1,:);
    
    for inh = 1:length(inhOffset)
        
        offset = inhOffset(inh);

            iExcite = zeros(totalCompartments,miscParams.time);
            iExcite(spine,250:300) = stimLvl;
            iInhibit = zeros(totalCompartments,miscParams.time);
            iInhibit(shaft,(250+offset):(300+offset)) = inhRange;
            
            input.excitaiton = iExcite;
            input.inhibition = iInhibit;
            
            voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
            yesInhibition(exc,inh,:) = voltages(1,:);
            
            inhDiff(exc,inh) = max(noInhibition(exc,:)) - max(yesInhibition(exc,inh,:));
            inhDiffAUC(exc,inh) = sum(noInhibition(exc,:)+60) - sum(yesInhibition(exc,inh,:)+60);
 
    end
    display(['Repeat ' num2str(r) ', stim level ' num2str(stimLvl)]);
end

inhEffect(r,:,:) = inhDiff;
inhEffectAUC(r,:,:) = inhDiffAUC;

end
tend = toc(tstart)

%% Process data
meanEffect = squeeze(mean(inhEffect,1));
rowMax = repmat(max(meanEffect,[],2),1,length(inhOffset));
normEffect = meanEffect./rowMax;

meanEffectAUC = squeeze(mean(inhEffectAUC,1));
rowMax = repmat(max(meanEffectAUC,[],2),1,length(inhOffset));
normEffectAUC = meanEffectAUC./rowMax;

%% Plot data
figure;
imagesc(meanEffect);
xlabel('Inhibition offset (ms)');
xlist = 1:5:length(inhOffset);
xticks(xlist);
xticklabels(inhOffset(xlist)./10);
ylabel('Excitation level');
ylist = 10:20:length(stimRange);
yticks(ylist);
yticklabels(stimRange(ylist));
cb = colorbar;
ylabel(cb,'Effect of inhibition')
title('Effect of inhibition timing (max voltage)');

figure;
imagesc(meanEffectAUC);
xlabel('Inhibition offset (ms)');
xlist = 1:5:length(inhOffset);
xticks(xlist);
xticklabels(inhOffset(xlist)./10);
ylabel('Excitation level');
ylist = 10:20:length(stimRange);
yticks(ylist);
yticklabels(stimRange(ylist));
cb = colorbar;
ylabel(cb,'Effect of inhibition')
title('Effect of inhibition timing (area under curve)');