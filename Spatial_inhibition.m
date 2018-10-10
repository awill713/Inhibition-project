
clear;
close all;

%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Design simulation parameters

repeats = 5;
stimRange = 1:1:10;
inhibitionLevel = 5;

%% Choose locations of inhibition
%if you add new compartments to the list, update the function
%compartmentFinder below to tell it how to find that compartment within the
%dendritic tree

inhLocator.distalTuft = @(cID,connect) max(find(cID(2,:)==3)); %most distal apical tuft

inhLocator.middleApical = @(cID,connect) intersect(find(cID(2,:)==2),... %find the most distal apical shaft
    find(connect(min(find(cID(2,:)==3)),:)==1)); %that is connected to the most proximal apical tuft

inhLocator.proximalApical = @(cID,connect) min(find(cID(2,:)==2)); %most proximal apical shaft
inhLocator.soma = @(cID,connect) 1;
inhLocator.proximalBasal = @(cID,connect) min(find(cID(2,:)==1)); %most proximal basal shaft
inhLocator.distalBasal = @(cID,connect) max(find(cID(2,:)==1)); %most distal basal shaft

inhibitionLocations = fieldnames(inhLocator);

%% Choose locations of excitation
%if you add new compartments to the list, update the function
%compartmentFinder below to tell it how to find that compartment within the
%dendritic tree
excLocator.distalTuft = @(cID,connect) max(find(connect(max(find(cID(2,:)==3)),:))); %spine connected to the most distal apical tuft
excLocator.middleApical = @(cID,connect) max(find(connect(intersect(find(cID(2,:)==2),... %spine connected to the most distal apical shaft
    find(connect(min(find(cID(2,:)==3)),:)==1)),:))); %that is connected to the most proximal apical tuft

excLocator.proximalApical = @(cID,connect) max(find(connect(min(find(cID(2,:)==2)),:))); %spine connected to the most proximal apical shaft
excLocator.soma = @(cID,connect) 1;
excLocator.proximalBasal = @(cID,connect) max(find(connect(min(find(cID(2,:)==1)),:))); %spine connected to the most proximal basal shaft
excLocator.distalBasal = @(cID,connect) max(find(connect(max(find(cID(2,:)==1)),:))); %spine connected to the most distal basal shaft

excitationLocations = fieldnames(excLocator);

%% Run simulations (repeats) for each combination of excitation/inhibition

output = zeros(length(excitationLocations),length(inhibitionLocations),miscParams.time,repeats);
for r = 1:repeats
    r
    % Build dendritic arbor
    [connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
    totalCompartments = size(compartmentIDs,2);
    
    %find the various compartments to be excited/inhibited in various
    %combinations
    [excComp, inhComp] = compartmentFinder(excLocator, inhLocator, compartmentIDs,connectome);
    
    %matrix holding soma voltage trace for entire experiment for each
    %combination of excitation and inhibition location. To be averaged
    %after all repetitions.
    tempOutput = zeros(length(excComp), length(inhComp),4);
    
    for exc = 1:length(excComp)
        
        %compartment to receive excitation this round
        excChoice = excComp(exc);
        
        for inh = 1:length(inhComp)
            
            %compartment to receive inhibition this round
            inhChoice = inhComp(inh);
            
            maxV = zeros(1,length(stimRange));
            somaAUC = zeros(1,length(stimRange));
            inhV = zeros(1,length(stimRange));
            inhAUC = zeros(1,length(stimRange));
            
            %First run simulations of increasing excitatory intensity
            %without inhibition, to get the "baseline" nonlinearity index
            for s = 1:length(stimRange)
                
                stimLevel = stimRange(s);
                
                %implement excited/inhibited compartments in this run's
                %simulation design
                iExcite = zeros(totalCompartments,miscParams.time);
                iExcite(excChoice,150:200) = stimLevel;
                
                iInhibit = zeros(totalCompartments,miscParams.time);
                iInhibit(inhChoice,150:200) = 0;
                
                input.excitation = iExcite;
                input.inhibition = iInhibit;
                
                %run simulation
                voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
                
                maxV(s) = max(voltages(1,:));
                somaAUC(s) = sum(voltages(1,150:250)-condParams.vRest);
                
                
                iInhibit(inhChoice,150:200) = inhibitionLevel;
                input.inhibition = iInhibit;
                
                %run simulation
                voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
                
                inhV(s) = max(voltages(1,:));
                inhAUC(s) = sum(voltages(1,150:250)-condParams.vRest);
                
            end
            
            predictedMaxV = -60+(maxV(1,1)+60).*(stimRange);
            predictedAUC = somaAUC(1,1)*stimRange;
            
            AUCpredictedMaxV = sum(predictedMaxV - predictedMaxV(1));
            AUCpredictedAUC = sum(predictedAUC - predictedAUC(1));
            
            AUCbaselineMaxV = sum(maxV - predictedMaxV(1));
            AUCbaselineAUC = sum(somaAUC - predictedAUC(1));
            
            AUCinhibitionMaxV = sum(inhV - predictedMaxV(1));
            AUCinhibitionAUC = sum(inhAUC - predictedAUC(1));
            
            baselineIndexMaxV  = AUCbaselineMaxV / AUCpredictedMaxV;
            baselineIndexAUC = AUCbaselineAUC / AUCpredictedAUC;
            inhibitionIndexMaxV = AUCinhibitionMaxV / AUCpredictedMaxV;
            inhibitionIndexAUC = AUCinhibitionAUC / AUCpredictedAUC;
            
            tempOutput(exc,inh,1) = baselineIndexMaxV;
            tempOutput(exc,inh,2) = baselineIndexAUC;
            tempOutput(exc,inh,3) = inhibitionIndexMaxV;
            tempOutput(exc,inh,4) = inhibitionIndexAUC;
            
        end
    end
end


%%

for i = 1:size(tempOutput,1)
    figure;
    for j = 1:size(tempOutput,2)
        subplot(1,size(tempOutput,2),j);hold on;
        bar(1,tempOutput(i,j,1));
        bar(2,tempOutput(i,j,3));
        xlabel(['Inhibition ' inhibitionLocations{j}]);
        if j == 1
            ylabel(['Excitation ' excitationLocations{i}]);
        end
    
    end
end

%{

%%

%store tempOutput (somatic voltage trace for each combination of
%excitation/inhibition) in ultimate variable "output." Could choose to
%directly store it in "output" without tempOutput, but thought this was
%a bit cleaner?
% output(:,:,:,r) = tempOutput;


somaTraceMean = mean(output,4);
somaTraceSTD = std(output,[],4);

%% Plot results

figureTime = 100:400;
for fExc = 1:size(somaTraceMean,1)
    figure(fExc)
%     figure(2*fExc-1);  hold on;
    subplot(2,2,1); hold on;
    for fInh = 1:size(somaTraceMean,2)
        plot(figureTime,squeeze(somaTraceMean(fExc,fInh,figureTime)));
    end
    xlabel('Time');
    ylabel('Soma voltage (mv)');
    title(['Excitation onto ' excitationLocations(fExc)]);
    legend(excitationLocations);
    
%     figure(2*fExc);
    subplot(2,2,2); hold on;
    traceMat = squeeze(somaTraceMean(fExc,:,figureTime));
    maxVector = max(squeeze(somaTraceMean(fExc,:,figureTime)),[],2);
    imagesc((traceMat+60)./(maxVector+60));
    xlabel('Time');
    xticks(1:100:size(traceMat,2));
    xticklabels(figureTime(1):100:figureTime(end))
    yticks(1:size(somaTraceMean,2));
    yticklabels(inhibitionLocations);
    title(['Excitation onto ' excitationLocations(fExc)]);
    colorbar
    
    subplot(2,2,3); hold on;
    [vMax maxInd] = max(squeeze(somaTraceMean(fExc,:,:)),[],2);
    plot(1:size(somaTraceMean,2),vMax);
    vMin = min(squeeze(somaTraceMean(fExc,:,:)),[],2);
    plot(1:size(somaTraceMean,2),vMin);
    xticks(1:size(somaTraceMean,2));
    xticklabels(inhibitionLocations);
    ylabel('Soma voltage');
    legend('Max voltage','Min voltage');
    title('Max and min soma voltage')
    
    subplot(2,2,4);
    plot(1:size(somaTraceMean,2),maxInd-150);
    xticks(1:size(somaTraceMean,2));
    xticklabels(inhibitionLocations);
    ylabel('Time relative to excitation impulse');
    title('Latency to peak voltage');
    
    %also plot duration and integration of deviation
end

ultimately though, after verifying that these inhibition results are
satisfactory and well understood, we want to compare the effects of
inhibition location on the sub/supralinear summation properties, which
entails cycling through a range of excitation intensities

%}

%%
    
function [excOut, inhOut] = compartmentFinder(excIn, inhIn, cID, connect)
    excOut = zeros(1,length(fieldnames(excIn)));
    inhOut = zeros(1,length(fieldnames(inhIn)));
    
    eFunc = fieldnames(excIn);
    for e = 1:length(eFunc)
        excOut(e) = feval(excIn.(eFunc{e}),cID,connect);
    end
    
    iFunc = fieldnames(inhIn);
    for i = 1:length(iFunc)
        inhOut(i) = feval(inhIn.(iFunc{i}),cID,connect);
    end
    
end
