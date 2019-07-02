
%%CHANGE INHIBITION LEVEL, explore parameter space, talk with Maria

clear;
% close all;

%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Design simulation parameters

repeats = 100;
stimRange = 1:1:30;
inhibitionLevel = [5 30 100 150]; %range of 20-40 seems like a working range?

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
inhCount = length(inhibitionLocations);

%% Choose locations of excitation
%if you add new compartments to the list, update the function
%compartmentFinder below to tell it how to find that compartment within the
%dendritic tree
excLocator.distalTuft = @(cID,connect) max(find(connect(max(find(cID(2,:)==3)),:))); %spine connected to the most distal apical tuft
excLocator.middleApical = @(cID,connect) max(find(connect(intersect(find(cID(2,:)==2),... %spine connected to the most distal apical shaft
    find(connect(min(find(cID(2,:)==3)),:)==1)),:))); %that is connected to the most proximal apical tuft

excLocator.proximalApical = @(cID,connect) max(find(connect(min(find(cID(2,:)==2)),:))); %spine connected to the most proximal apical shaft
% excLocator.soma = @(cID,connect) 1;
excLocator.proximalBasal = @(cID,connect) max(find(connect(min(find(cID(2,:)==1)),:))); %spine connected to the most proximal basal shaft
excLocator.distalBasal = @(cID,connect) max(find(connect(max(find(cID(2,:)==1)),:))); %spine connected to the most distal basal shaft

excitationLocations = fieldnames(excLocator);
excCount = length(excitationLocations);

%% Run simulations (repeats) for each combination of excitation/inhibition/stimlevel

withoutInhibition = zeros(repeats,length(excitationLocations),length(inhibitionLocations),length(stimRange),length(inhibitionLevel),miscParams.time);
withInhibition = zeros(repeats,length(excitationLocations),length(inhibitionLocations),length(stimRange),length(inhibitionLevel),miscParams.time);
tstart = tic;
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
                
                for inhLvl = 1:length(inhibitionLevel)
                    iLevel = inhibitionLevel(inhLvl);
                    
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
                    
                    withoutInhibition(r,exc,inh,s,inhLvl,:) = voltages(1,:);
                    %                 maxV(s) = max(voltages(1,:));
                    %                 somaAUC(s) = sum(voltages(1,150:250)-condParams.vRest);
                    
                    
                    iInhibit(inhChoice,150:200) = iLevel;
                    input.inhibition = iInhibit;
                    
                    %run simulation
                    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
                    
                    withInhibition(r,exc,inh,s,inhLvl,:) = voltages(1,:);
                    %                 inhV(s) = max(voltages(1,:));
                    %                 inhAUC(s) = sum(voltages(1,150:250)-condParams.vRest);
                end
            end
            
            %             predictedMaxV = -60+(maxV(1,1)+60).*(stimRange);
            %             predictedAUC = somaAUC(1,1)*stimRange;
            %
            %             AUCpredictedMaxV = sum(predictedMaxV - predictedMaxV(1));
            %             AUCpredictedAUC = sum(predictedAUC - predictedAUC(1));
            %
            %             AUCbaselineMaxV = sum(maxV - predictedMaxV(1));
            %             AUCbaselineAUC = sum(somaAUC - predictedAUC(1));
            %
            %             AUCinhibitionMaxV = sum(inhV - predictedMaxV(1));
            %             AUCinhibitionAUC = sum(inhAUC - predictedAUC(1));
            %
            %             baselineIndexMaxV  = AUCbaselineMaxV / AUCpredictedMaxV;
            %             baselineIndexAUC = AUCbaselineAUC / AUCpredictedAUC;
            %             inhibitionIndexMaxV = AUCinhibitionMaxV / AUCpredictedMaxV;
            %             inhibitionIndexAUC = AUCinhibitionAUC / AUCpredictedAUC;
            %
            %             tempOutput(exc,inh,1) = baselineIndexMaxV;
            %             tempOutput(exc,inh,2) = baselineIndexAUC;
            %             tempOutput(exc,inh,3) = inhibitionIndexMaxV;
            %             tempOutput(exc,inh,4) = inhibitionIndexAUC;
            
        end
    end
    telapsed = toc(tstart)
end

%% Analyze data
% color = {'k','b','r','g','m','c';'k--','b--','r--','g--','m--','c--'}
% figure; hold on;
% for i = 6: 6
%     plot(cumsum(squeeze(withoutInhibition(5,5,i,7,:))+60),'k');
%     plot(cumsum(squeeze(withInhibition(5,5,i,7,:))+60),'k--');
% end
% 
% % color = {'k','b','r','g','m','c';'k--','b--','r--','g--','m--','c--'}
% color = {'k','b','b','k','m','c';'k--','b--','b--','k--','b--','c--'}
% 
% figure; hold on;
% for i = 6:6
%     plot(squeeze(mean(withoutInhibition(:,3,i,7,:))),'k');
%     plot(squeeze(mean(withInhibition(:,3,i,7,:))),'k--');
% end




% for ii = 1:length(inhibitionLevel)

%exc 7 inh 30, exc 13/8/7 inh 100, exc 13/20/21 inh 150
for ii = 1:1
    inhIntense = inhibitionLevel(ii);
% intensitiesToPlot = [1 4 5 6 7 8 10];
intensitiesToPlot = [1 5 6 7 8 11 12 13 19 20 21 22];
% intensitiesToPlot = 1:1:30;
for i = 1:length(intensitiesToPlot)
    intense = intensitiesToPlot(i);
    figure;
%     noInh = squeeze(mean(max(withoutInhibition(:,:,:,intense,:),[],5),1));
%     yesInh = squeeze(mean(max(withInhibition(:,:,:,intense,:),[],5),1));
    
    noInh = squeeze(max(mean(withoutInhibition(:,:,:,intense,ii,:),1),[],6));
    yesInh = squeeze(max(mean(withInhibition(:,:,:,intense,ii,:),1),[],6));
    
%     noInh = squeeze(max(cumsum(mean(withoutInhibition(:,:,:,intense,ii,:),1)-condParams.vRest,6),[],6));
%     yesInh = squeeze(max(cumsum(mean(withInhibition(:,:,:,intense,ii,:),1)-condParams.vRest,6),[],6));
    
    diff = noInh - yesInh;
    rowMax = repmat(max(diff,[],2),1,inhCount);
    imagesc((noInh - yesInh)./rowMax);
%     imagesc(noInh - yesInh);
    c = colorbar;
    xlabel('Inhibition location');
    xticks(1:1:inhCount);
    xticklabels(inhibitionLocations);
    ylabel('Excitation location');
    yticks(1:1:excCount);
    yticklabels(excitationLocations);
    ylabel(c,'Normalized difference in max voltage with inhibition')
    title(['Exc = ' num2str(intense) ', inh = ' num2str(inhIntense)]);
end
end

% indices = zeros(excCount,inhCount);
% indicesInh = zeros(excCount,inhCount);
% traject = zeros(excCount,inhCount,repeats,2,length(stimRange));
% rSquared = zeros(excCount,inhCount,2);
% for e = 1:excCount
%     
%     for i = 1:inhCount
%         
%     tempIndex = zeros(2,repeats);
%     
%     for r = 1:repeats
%         predicted = (-60 + (max(withoutInhibition(r,e,i,1,:))+60) * stimRange);
%         AUCpredicted = sum(predicted - condParams.vRest);
%         
%         traject(e,i,r,1,:) = squeeze(max(withoutInhibition(r,e,i,:,:),[],5));
%         AUCnoInh = sum(squeeze(max(withoutInhibition(r,e,i,:,:),[],5)) - condParams.vRest);
%         nonLinIndex = (AUCnoInh - AUCpredicted) / AUCpredicted;
%         
%         traject(e,i,r,2,:) = squeeze(max(withInhibition(r,e,i,:,:),[],5));
%         AUCinh = sum(squeeze(max(withInhibition(r,e,i,:,:),[],5)) - condParams.vRest);
%         nonLinIndexInh = (AUCinh - AUCpredicted) / AUCpredicted;
%         
%         tempIndex(1,r) = nonLinIndex;
%         tempIndex(2,r) = nonLinIndexInh;
%     end
%     
%     meanTraject = squeeze(mean(traject,3));
%     
%     pf = polyfit(stimRange,squeeze(meanTraject(e,i,1,:))',1);
%     pv = polyval(pf,stimRange);
%     Bbar = mean(squeeze(meanTraject(e,i,1,:)));
%     SStot = sum((squeeze(meanTraject(e,i,1,:))' - Bbar).^2);
%     SSres = sum((squeeze(meanTraject(e,i,1,:))' - pv).^2);
%     R2 = 1 - SSres/SStot;
%     rSquared(e,i,1) = R2;
%     
%     pf = polyfit(stimRange,squeeze(meanTraject(e,i,2,:))',1);
%     pv = polyval(pf,stimRange);
%     Bbar = mean(squeeze(meanTraject(e,i,2,:)));
%     SStot = sum((squeeze(meanTraject(e,i,2,:))' - Bbar).^2);
%     SSres = sum((squeeze(meanTraject(e,i,2,:))' - pv).^2);
%     R2 = 1 - SSres/SStot;
%     rSquared(e,i,2) = R2;
%     
%     indices(e,i) = mean(tempIndex(1,:));
%     indicesInh(e,i) = mean(tempIndex(2,:));
%     end
% end
% figure;colormap jet;
% imagesc(indices);
% colorbar;
% title('Nonlinearity index');
% xticks(1:1:inhCount);
% xticklabels(inhibitionLocations);
% ylabel('Excitation location');
% yticks(1:1:excCount);
% yticklabels(excitationLocations);
% 
% figure;colormap jet;
% imagesc(indicesInh - indices);
% c = colorbar;
% title('Difference in nonlinearity index');
% xlabel('Inhibition location');
% xticks(1:1:inhCount);
% xticklabels(inhibitionLocations);
% ylabel('Excitation location');
% yticks(1:1:excCount);
% yticklabels(excitationLocations);
% ylabel(c,'Nonlinearity index with inhibition - without inhibition');
% 
% figure;colormap jet;
% imagesc(indicesInh ./ indices);
% c = colorbar;
% title('Ratio of nonlinearity index');
% xlabel('Inhibition location');
% xticks(1:1:inhCount);
% xticklabels(inhibitionLocations);
% ylabel('Excitation location');
% yticks(1:1:excCount);
% yticklabels(excitationLocations);
% ylabel(c,'Nonlinearity index with inhibition / without inhibition');
% 
% figure;colormap jet;
% imagesc(squeeze(rSquared(:,:,1)));
% c = colorbar;
% title('R-squared without inhibition)');
% xlabel('Inhibition location');
% xticks(1:1:inhCount);
% xticklabels(inhibitionLocations);
% ylabel('Excitation location');
% yticks(1:1:excCount);
% yticklabels(excitationLocations);
% ylabel(c,'R-squared value');
% 
% figure;colormap jet;
% imagesc(squeeze(rSquared(:,:,2)) - squeeze(rSquared(:,:,1)));
% c = colorbar;
% title('Difference in R-squared (Inh - noInh)');
% xlabel('Inhibition location');
% xticks(1:1:inhCount);
% xticklabels(inhibitionLocations);
% ylabel('Excitation location');
% yticks(1:1:excCount);
% yticklabels(excitationLocations);
% ylabel(c,'R-squared difference (with inhibition - without inhibition)');


%%

% for i = 1:size(tempOutput,1)
%     figure;
%     for j = 1:size(tempOutput,2)
%         subplot(1,size(tempOutput,2),j);hold on;
%         bar(1,tempOutput(i,j,1));
%         bar(2,tempOutput(i,j,3));
%         xlabel(['Inhibition ' inhibitionLocations{j}]);
%         if j == 1
%             ylabel(['Excitation ' excitationLocations{i}]);
%         end
%     
%     end
% end

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
