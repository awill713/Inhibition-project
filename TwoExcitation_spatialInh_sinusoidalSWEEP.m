
clear;
% close all;

stimMeanRange = 1:0.1:10;
inhRange = 0:10:40;
repeats = 10;

outputNorm = zeros(length(stimMeanRange),length(inhRange));
outputStand = zeros(length(stimMeanRange),length(inhRange));

%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Design simulation parameters
for stim = 1:length(stimMeanRange)
    stimValue = stimMeanRange(stim);
    
    for inh = 1:length(inhRange)
        inhValue = inhRange(inh);
        
        excitationMean = stimValue; %originally 15
        excitationAmplitude = 2*stimValue; %originally 10
        frequency1 = 0.01; %actual frequency in real time is this/2pi (multiply by 10000 for Hz)
        frequency2 = 0.006; %was 0.03
        
        inhibitionIntensity = inhValue;
        
        %% Storage variables
        spineList = zeros(repeats,2);
        shaftList = zeros(repeats,2);
        branchCompartment = zeros(1,repeats);
        
        excPower = zeros(repeats,3,50001);
        inhPower = zeros(repeats,2,50001);
        powerMaxIndices = zeros(repeats,2);
        powerRatio = zeros(repeats,5);
        normRatio = zeros(repeats,5);
        standRatio = zeros(repeats,5);
        
        %% Run simulations
        
        for r = 1:repeats
            display(['Stim = ' num2str(stimValue) ', inh = ' num2str(inhValue) ', repeat = ' num2str(r)]);
            
            % Build arbor
            [connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
            g = graph(connectome);
            graphMaster(r).arbor = g;
            totalCompartments = length(compartmentIDs);
            shaftCmpt = compartmentIDs(1,find(compartmentIDs(2,:)~=4));
            spineCmpt = compartmentIDs(1,find(compartmentIDs(2,:)==4));
            
            spine1 = datasample(spineCmpt,1);
            badList = [1 spine1]; %spine2 can't be soma or spine1
            spine2 = datasample(setdiff(spineCmpt,badList),1);
            path1 = shortestpath(g,spine1,1);path1 = path1(2:end);
            path2 = shortestpath(g,spine2,1);path2 = path2(2:end);
            while min(ismember(path1,path2))==1 || min(ismember(path2,path1))==1 %|| max(intersect(path1,path2))==1
                display(['Repicking spines']);
                
                spine1 = datasample(spineCmpt,1);
                badList = [1 spine1];
                spine2 = datasample(setdiff(spineCmpt,badList),1);
                
                path1 = shortestpath(g,spine1,1);path1 = path1(2:end);
                path2 = shortestpath(g,spine2,1);path2 = path2(2:end);
            end
            spineList(r,:) = [spine1 spine2];
            shaft1 = find(connectome(spine1,:)==1);
            shaft2 = find(connectome(spine2,:)==1);
            shaftList(r,:) = [shaft1 shaft2];
            % branchCompartment(1,r) = max(intersect(path1,path2));
            % pathCompartments{r,1} = path1;
            % pathCompartments{r,2} = path2;
            % path1Length(r) = length(path1);
            % path2Length(r) = length(path2);
            
            %% Design stimulus
            spine1Stimulus = zeros(1,miscParams.time);
            spine1Stimulus(1,:) = excitationAmplitude/2*cos((1:miscParams.time)*frequency1)+excitationMean;
            spine2Stimulus = zeros(1,miscParams.time);
            spine2Stimulus(1,:) = excitationAmplitude/2*cos((1:miscParams.time)*frequency2)+excitationMean;
            
            %% Define pure excitatory dynamics
            iExcite = zeros(totalCompartments,miscParams.time);
            iExcite(spine1,:) = spine1Stimulus;
            iInhibit = zeros(totalCompartments,miscParams.time);
            input.excitation = iExcite;
            input.inhibition = iInhibit;
            voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
            % figure;hold on;
            % plot(voltages(1,:));
            % % plot(spine1Stimulus);
            % yyaxis right;
            % plot(voltages(spine1,:));
            % title('Spine 1 excitation');
            v1 = voltages(1,:);
            [excPowerTemp(1,1:50001) freq] = periodogram(v1,hamming(length(v1)),100000,10000);
            temp = abs(freq - (frequency1*10000/(2*pi)));freq1Index = find(temp==min(temp));
            spine1Max = find(excPowerTemp(1,:)==max(excPowerTemp(1,freq1Index-20:freq1Index)));
            
            iExcite = zeros(totalCompartments,miscParams.time);
            iExcite(spine2,:) = spine2Stimulus;
            iInhibit = zeros(totalCompartments,miscParams.time);
            input.excitation = iExcite;
            input.inhibition = iInhibit;
            voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
            % figure;hold on;
            % plot(voltages(1,:));
            % % plot(spine2Stimulus);
            % yyaxis right;
            % plot(voltages(spine2,:));
            % title('Spine 2 excitation');
            v2 = voltages(1,:);
            excPowerTemp(2,1:50001) = periodogram(v2,hamming(length(v2)),100000,10000);
            temp = abs(freq - (frequency2*10000/(2*pi)));freq2Index = find(temp==min(temp));
            spine2Max = find(excPowerTemp(2,:)==max(excPowerTemp(2,freq2Index-20:freq2Index)));
            
            iExcite = zeros(totalCompartments,miscParams.time);
            iExcite(spine1,:) = spine1Stimulus;
            iExcite(spine2,:) = spine2Stimulus;
            iInhibit = zeros(totalCompartments,miscParams.time);
            input.excitation = iExcite;
            input.inhibition = iInhibit;
            voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
            % figure;hold on;
            % plot(voltages(1,:));
            % yyaxis right;
            % plot(voltages(spine1,:));
            % plot(voltages(spine2,:));
            % title('Spine 1 and 2 excitation');
            v12 = voltages(1,:);
            excPowerTemp(3,1:50001) = periodogram(v12,hamming(length(v12)),100000,10000);
            
            %
            % figure;hold on;
            % plot(excPower(1,:));
            % plot(excPower(2,:));
            % plot(excPower(3,:));
            % line([spine1Max spine1Max],[0 1],'Color','black','LineStyle','--');
            % line([spine2Max spine2Max],[0 1],'Color','black','LineStyle','--');
            % xlim([50 200]);
            % legend('Spine 1','Spine 2','Both spines');
            % title('Excitatory dynamics');
            
            % figure;hold on;
            % plot(v1); plot(v2); plot(v12);
            % legend({'V1','V2','V12'});
            
            %% Sweep through shafts to inhibit
            
            inhPowerTemp = zeros(length(shaftCmpt),50001);
            
            % for s = 1:length(shaftCmpt)
            shaftsToInhibit = [shaft1 shaft2];
            for sha = 1:length(shaftsToInhibit);
                s = shaftsToInhibit(sha);
                shaft = shaftCmpt(s);
                
                %Excite spines 1 and 2
                iExcite = zeros(totalCompartments,miscParams.time);
                iExcite(spine1,:) = spine1Stimulus;
                iExcite(spine2,:) = spine2Stimulus;
                iInhibit = zeros(totalCompartments,miscParams.time);
                iInhibit(shaft,:) = inhibitionIntensity;
                input.excitation = iExcite;
                input.inhibition = iInhibit;
                voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
                
                inhPowerTemp(s,:) = periodogram(voltages(1,:),hamming(length(voltages(1,:))),100000,10000);
            end
            
            %% Process data
            pow(1) = excPowerTemp(1,spine1Max) / excPowerTemp(1,spine2Max);
            pow(2) = excPowerTemp(2,spine1Max) / excPowerTemp(2,spine2Max);
            pow(3) = excPowerTemp(3,spine1Max) / excPowerTemp(3,spine2Max);
            pow(4) = inhPowerTemp(shaft1,spine1Max) / inhPowerTemp(shaft1,spine2Max);
            pow(5) = inhPowerTemp(shaft2,spine1Max) / inhPowerTemp(shaft2,spine2Max);
            normPow = pow./pow(3);
            standPow = (2*pow-pow(1)-pow(2))./(pow(1)+pow(2));
            
            %% Store data
            excPower(r,:,:) = [excPowerTemp(1,:); excPowerTemp(2,:); excPowerTemp(3,:)];
            inhPower(r,:,:) = inhPowerTemp([shaft1 shaft2],:);
            powerMaxIndices(r,:) = [spine1Max spine2Max];
            powerRatio(r,:) = pow;
            normRatio(r,:) = normPow;
            standRatio(r,:) = standPow;
            
        end
        
        %% Average across repeats
        meanExcPower = squeeze(mean(excPower));
        meanInhPower = squeeze(mean(inhPower));
        meanPowerRatio = mean(powerRatio);
        meanNormRatio = mean(normRatio);
        meanStandRatio = mean(standRatio);
        meanIndices = mean(powerMaxIndices);
        
        %%
        
        outputNorm(stim,inh) = meanNormRatio(1,5) - meanNormRatio(1,4);
        outputStand(stim,inh) = meanStandRatio(1,5) - meanStandRatio(1,4);
    end
end

figure;
imagesc(outputNorm);
c = colorbar;
xlabel('Inhibition level');
xticks([1 length(inhRange)]);
xticklabels([inhRange(1), inhRange(end)]);
ylabel('Excitation level');
yticks([1 length(stimMeanRange)]);
yticklabels([stimMeanRange(1), stimMeanRange(end)]);
ylabel(c,'Difference in freq1:freq2 (inh2 - inh1)')
title('Normalized around 1');

figure;
imagesc(outputStand);
c = colorbar;
xlabel('Inhibition level');
xticks([1 length(inhRange)]);
xticklabels([inhRange(1), inhRange(end)]);
ylabel('Excitation level');
yticks([1 length(stimMeanRange)]);
yticklabels([stimMeanRange(1), stimMeanRange(end)]);
ylabel(c,'Difference in freq1:freq2 (inh2 - inh1)')
title('Standardized (-1 to 1)')
