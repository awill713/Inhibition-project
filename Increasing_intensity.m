%This entire script simulates increasing excitatory input onto a single
%shaft (apical or basal) or spine (apical or basal)

clear;
% close all;

%% Load parameters
[miscParams dendParams condParams ] = loadParameters;

%% Build dendritic arbor
[connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
totalCompartments = size(compartmentIDs,2);

%%
% THIS SECTION COMPARES INCREASING EXCITATION INTENSITY, COMPARING
% SUB/SUPRALINEAR SUMMATION OF SPINES AND SHAFTS IN APICAL OR BASAL
% DENDRITES

stimOn = 150;
stimOff = 200;
maxIntensity = 40;

%% Apical shaft - specify excitatory and inhibitory input
iExcite = zeros(maxIntensity,totalCompartments,miscParams.time);
iInhibit = zeros(maxIntensity,totalCompartments,miscParams.time);

stimCompartment = max(find(compartmentIDs(2,:)==3)); %most distal shaft
for intense = 1:maxIntensity
    iExcite(intense,stimCompartment,stimOn:stimOff) = intense;
end

maxV = zeros(1,size(iExcite,1));
somaAUC = zeros(1,size(iExcite,1));

%% Apical shaft - run model
for exp = 1:size(iExcite,1)
    exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(exp) = max(voltages(1,:));
    somaAUC(exp) = sum(voltages(1,150:250)-condParams.vRest);
end

AUCpredicted = sum((-60+(maxV(1,1)+60).*(1:1:maxIntensity)) - condParams.vRest);
AUCmaxV = (sum(maxV-condParams.vRest) - AUCpredicted) / AUCpredicted;
f1 = figure;plot(maxV);
hold on;
plot(-60+(maxV(1,1)+60).*(1:1:maxIntensity));
legend('Observed','Predicted');
title(['Apical shaft excitation is (sub)linear, nonlinearity index = ' num2str(AUCmaxV)]);
xlabel('Excitation intensity (nS)');
ylabel('Maximum soma membrane potential (mV)');

AUCpredicted = sum(somaAUC(1,1).*(1:1:maxIntensity));
AUCsomaAUC = (sum(somaAUC) - AUCpredicted) / AUCpredicted;
f2 = figure;plot(somaAUC);
hold on;
plot(somaAUC(1,1).*(1:1:maxIntensity));
legend('Observed','Predicted');
title(['Apical shaft area under curve, nonlinearity index = ' num2str(AUCsomaAUC)]);
xlabel('Excitation intensity (nS)');
ylabel('Area under curve (mV * s)');

iExcite = zeros(maxIntensity,totalCompartments,miscParams.time);
iInhibit = zeros(maxIntensity,totalCompartments,miscParams.time);

%% Apical spine - specify excitatory and inhibitory input
stimCompartment = max(find(compartmentIDs(2,:)==4)); %most distal apical spine
for intense = 1:maxIntensity
    iExcite(intense,stimCompartment,stimOn:stimOff) = intense;
end

maxV = zeros(1,size(iExcite,1));
somaAUC = zeros(1,size(iExcite,1));

%% Apical spine - run model
for exp = 1:size(iExcite,1)
    exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(exp) = max(voltages(1,:));
    somaAUC(exp) = sum(voltages(1,150:250)-condParams.vRest);
end

AUCpredicted = sum((-60+(maxV(1,1)+60).*(1:1:maxIntensity)) - condParams.vRest);
AUCmaxV = (sum(maxV-condParams.vRest) - AUCpredicted) / AUCpredicted;
f3 = figure;plot(maxV);
hold on;
plot(-60+(maxV(1,1)+60).*(1:1:maxIntensity));
legend('Observed','Predicted');
title(['Apical spine excitation is supralinear, nonlinearity index = ' num2str(AUCmaxV)]);
xlabel('Excitation intensity (nS)');
ylabel('Maximum soma membrane potential (mV)');

AUCpredicted = sum(somaAUC(1,1).*(1:1:maxIntensity));
AUCsomaAUC = (sum(somaAUC) - AUCpredicted) / AUCpredicted;
f4 = figure;plot(somaAUC);
hold on;
plot(somaAUC(1,1).*(1:1:maxIntensity));
legend('Observed','Predicted');
title(['Apical spine area under curve, nonlinearity index = ' num2str(AUCsomaAUC)]);
xlabel('Excitation intensity (nS)');
ylabel('Area under curve (mV * s)');

%% Basal shaft - specify excitatory and inhibitory input
iExcite = zeros(maxIntensity,totalCompartments,miscParams.time);
iInhibit = zeros(maxIntensity,totalCompartments,miscParams.time);

stimCompartment = max(find(compartmentIDs(2,:)==1)); %most distal basal shaft
for intense = 1:maxIntensity
    iExcite(intense,stimCompartment,stimOn:stimOff) = intense;
end

maxV = zeros(1,size(iExcite,1));
somaAUC = zeros(1,size(iExcite,1));

%% Basal shaft - run model
for exp = 1:size(iExcite,1)
    exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(exp) = max(voltages(1,:));
    somaAUC(exp) = sum(voltages(1,150:250)-condParams.vRest);
end

AUCpredicted = sum((-60+(maxV(1,1)+60).*(1:1:maxIntensity)) - condParams.vRest);
AUCmaxV = (sum(maxV-condParams.vRest) - AUCpredicted) / AUCpredicted;
f5 = figure;plot(maxV);
hold on;
plot(-60+(maxV(1,1)+60).*(1:1:maxIntensity));
legend('Observed','Predicted');
title(['Basal shaft excitation is (sub)linear, nonlinearity index = ' num2str(AUCmaxV)]);
xlabel('Excitation intensity (nS)');
ylabel('Maximum soma membrane potential (mV)');

AUCpredicted = sum(somaAUC(1,1).*(1:1:maxIntensity));
AUCsomaAUC = (sum(somaAUC) - AUCpredicted) / AUCpredicted;
f6 = figure;plot(somaAUC);
hold on;
plot(somaAUC(1,1).*(1:1:maxIntensity));
legend('Observed','Predicted');
title(['Basal shaft area under curve, nonlinearity index = ' num2str(AUCsomaAUC)]);
xlabel('Excitation intensity (nS)');
ylabel('Area under curve (mV * s)');

%% Basal spine - specify excitatory and inhibitory input
iExcite = zeros(maxIntensity,totalCompartments,miscParams.time);
iInhibit = zeros(maxIntensity,totalCompartments,miscParams.time);

stimCompartment = max(find(connectome(max(find(compartmentIDs(2,:)==1)),:))); %most distal basal spine
for intense = 1:maxIntensity
    iExcite(intense,stimCompartment,stimOn:stimOff) = intense;
end

maxV = zeros(1,size(iExcite,1));
somaAUC = zeros(1,size(iExcite,1));

%% Basal spine - run model
for exp = 1:size(iExcite,1)
    exp
    input.excitation = squeeze(iExcite(exp,:,:));
    input.inhibition = squeeze(iInhibit(exp,:,:));
    voltages = runSimulation(conductanceMat, compartmentIDs, input, condParams);
    maxV(exp) = max(voltages(1,:));
    somaAUC(exp) = sum(voltages(1,150:250)-condParams.vRest);
end

AUCpredicted = sum((-60+(maxV(1,1)+60).*(1:1:maxIntensity)) - condParams.vRest);
AUCmaxV = (sum(maxV-condParams.vRest) - AUCpredicted) / AUCpredicted;
f7 = figure;plot(maxV);
hold on;
plot(-60+(maxV(1,1)+60).*(1:1:maxIntensity));
legend('Observed','Predicted');
title(['Basal spine excitation is supralinear, nonlinearity index = ' num2str(AUCmaxV)]);
xlabel('Excitation intensity (nS)');
ylabel('Maximum soma membrane potential (mV)');

AUCpredicted = sum(somaAUC(1,1).*(1:1:maxIntensity));
AUCsomaAUC = (sum(somaAUC) - AUCpredicted) / AUCpredicted;
f8 = figure;plot(somaAUC);
hold on;
plot(somaAUC(1,1).*(1:1:maxIntensity));
legend('Observed','Predicted');
title(['Basal spine area under curve, nonlinearity index = ' num2str(AUCsomaAUC)]);
xlabel('Excitation intensity (nS)');
ylabel('Area under curve (mV * s)');
