
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% MAKE misc.totalTime=500 IN loadParameters!!! %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This script looks at the ability of inhibition to filter out sinusoidally
%oscillating inputs on spine1 and spine2. It looks at the power of the two
%oscillatory frequencies felt in the soma, and compares it with and without
%inhibition.

%%
clear;
% close all;

repeats = 10;

%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Design simulation parameters

%Oscillatory excitatory input parameters
excitationMean = 10; %originally 15
excitationAmplitude = 10; %originally 10
frequency1 = 0.01; %actual frequency in real time is this/2pi (multiply by 10000 for Hz)
frequency2 = 0.006; %was 0.03

inhibitionIntensity = 40;

%% Storage variables
spineList = zeros(repeats,2);
shaftList = zeros(repeats,2);
branchCompartment = zeros(1,repeats);

excPower = zeros(repeats,3,50001); %3 for exciting spine1, spine2 and spine1/2. 50001 due to the resolution of periodogram function
inhPower = zeros(repeats,2,50001); %2 for inhibiting shaft1 and shaft2
powerMaxIndices = zeros(repeats,2); %power at frequency1 and frequency2
powerRatio = zeros(repeats,5); %exc 1, exc2, exc1/2, exc1/2 inh1, exc1/2 inh2
normRatio = zeros(repeats,5); %power ratio of exc1/2 = 1
standardizedRatio = zeros(repeats,5); %power ratio of exc1 = 1, power ratio of exc2 = -1

%% Run simulations

for r = 1:repeats
    r

%% Build arbor and choose compartments
[connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
g = graph(connectome);
graphMaster(r).arbor = g;
totalCompartments = length(compartmentIDs);
shaftCmpt = compartmentIDs(1,find(compartmentIDs(2,:)~=4));
spineCmpt = compartmentIDs(1,find(compartmentIDs(2,:)==4));

spine1 = datasample(spineCmpt,1);
badList = [spine1]; %spine2 can't be spine1
spine2 = datasample(setdiff(spineCmpt,badList),1);
path1 = shortestpath(g,spine1,1);path1 = path1(2:end);
path2 = shortestpath(g,spine2,1);path2 = path2(2:end);
while min(ismember(path1,path2))==1 || min(ismember(path2,path1))==1 %|| max(intersect(path1,path2))==1
    display(['Repicking spines']);
    
    spine1 = datasample(spineCmpt,1);
    badList = [spine1];
    spine2 = datasample(setdiff(spineCmpt,badList),1);

    path1 = shortestpath(g,spine1,1);path1 = path1(2:end);
    path2 = shortestpath(g,spine2,1);path2 = path2(2:end);
end
spineList(r,:) = [spine1 spine2];
shaft1 = find(connectome(spine1,:)==1);
shaft2 = find(connectome(spine2,:)==1);
shaftList(r,:) = [shaft1 shaft2];

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
v1 = voltages(1,:);
[excPowerTemp(1,1:50001) freq] = periodogram(v1,hamming(length(v1)),100000,10000);
temp = abs(freq - (frequency1*10000/(2*pi)));freq1Index = find(temp==min(temp));
spine1Max = find(excPowerTemp(1,:)==max(excPowerTemp(1,freq1Index-20:freq1Index))); %index in periodogram of frequency1

iExcite = zeros(totalCompartments,miscParams.time);
iExcite(spine2,:) = spine2Stimulus;
iInhibit = zeros(totalCompartments,miscParams.time);
input.excitation = iExcite;
input.inhibition = iInhibit;
voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
v2 = voltages(1,:);
excPowerTemp(2,1:50001) = periodogram(v2,hamming(length(v2)),100000,10000);
temp = abs(freq - (frequency2*10000/(2*pi)));freq2Index = find(temp==min(temp));
spine2Max = find(excPowerTemp(2,:)==max(excPowerTemp(2,freq2Index-20:freq2Index))); %index in periodogram of frequency2

iExcite = zeros(totalCompartments,miscParams.time);
iExcite(spine1,:) = spine1Stimulus;
iExcite(spine2,:) = spine2Stimulus;
iInhibit = zeros(totalCompartments,miscParams.time);
input.excitation = iExcite;
input.inhibition = iInhibit;
voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);

v12 = voltages(1,:);
excPowerTemp(3,1:50001) = periodogram(v12,hamming(length(v12)),100000,10000);

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

%ratio between power at freq1 and power at freq1
pow(1) = excPowerTemp(1,spine1Max) / excPowerTemp(1,spine2Max); %exc1
pow(2) = excPowerTemp(2,spine1Max) / excPowerTemp(2,spine2Max); %exc2
pow(3) = excPowerTemp(3,spine1Max) / excPowerTemp(3,spine2Max); %exc1/2
pow(4) = inhPowerTemp(shaft1,spine1Max) / inhPowerTemp(shaft1,spine2Max); %exc1/2 inh1
pow(5) = inhPowerTemp(shaft2,spine1Max) / inhPowerTemp(shaft2,spine2Max); %exc1/2 inh2

normPow = pow./pow(3); %power ratio of exc1/2 = 1
standPow = (2*pow-pow(1)-pow(2))./(pow(1)+pow(2)); %power ratio of exc1 = 1, ratio of exc2 = -1

%% Store data
excPower(r,:,:) = [excPowerTemp(1,:); excPowerTemp(2,:); excPowerTemp(3,:)];
inhPower(r,:,:) = inhPowerTemp([shaft1 shaft2],:);
powerMaxIndices(r,:) = [spine1Max spine2Max];
powerRatio(r,:) = pow;
normRatio(r,:) = normPow;
standardizedRatio(r,:) = standPow;

end

%% Average across repeats
meanExcPower = squeeze(mean(excPower));
meanInhPower = squeeze(mean(inhPower));
meanPowerRatio = mean(powerRatio);
meanNormRatio = mean(normRatio);
meanStandRatio = mean(standardizedRatio);
meanIndices = mean(powerMaxIndices);

%% Display results
figure;hold on;
plot(meanInhPower(1,:));
plot(meanInhPower(2,:));
plot(meanExcPower(3,:));
line([meanIndices(1) meanIndices(1)],[0 1],'Color','blue','LineStyle','--');
line([meanIndices(2) meanIndices(2)],[0 1],'Color','red','LineStyle','--');
xlim([51 201]);
xlist = 51:50:201;
xticks(xlist);
xticklabels(freq(xlist));
xlabel('Frequency (Hz');
ylabel('Power');
legend('Inh shaft 1','Inh shaft 2','No inhibition');
title('Excitation compared to inhibition');

figure;hold on;
plot(meanInhPower(1,:)-meanExcPower(3,:));
plot(meanInhPower(2,:)-meanExcPower(3,:));
line([meanIndices(1) meanIndices(1)],[-0.1 0.1],'Color','blue','LineStyle','--');
line([meanIndices(2) meanIndices(2)],[-0.1 0.1],'Color','red','LineStyle','--');
xlim([51 201]);
xlist = 51:50:201;
xticks(xlist);
xticklabels(freq(xlist));
xlabel('Frequency (Hz');
ylabel('\Delta power');
legend('Inhibit shaft1','Inhibit shaft2')
title('Difference in power with inhibition');

figure;bar([meanNormRatio(1,3); meanNormRatio(1,4); meanNormRatio(1,5)]);
% line([0.5 3.5],[meanNormRatio(1,1) meanNormRatio(1,1)],'LineStyle','--');
% line([0.5 3.5],[meanNormRatio(1,2) meanNormRatio(1,2)],'LineStyle','--');
xticklabels({'No inhibition','Inhibit shaft 1','Inhibit shaft 2'});
title('Signal 1 : signal 2 (no inh = 1)');

figure;bar([meanStandRatio(1,3); meanStandRatio(1,4); meanStandRatio(1,5)]);
line([0.5 3.5],[1 1],'LineStyle','--');
line([0.5 3.5],[-1 -1],'LineStyle','--');
xticklabels({'No inhibition','Inhibit shaft 1','Inhibit shaft 2'});
title('Signal 1 : signal 2 (-1 : 1)');