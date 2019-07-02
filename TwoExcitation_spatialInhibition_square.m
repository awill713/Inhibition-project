%This script looks at how the location of inhibition determines its ability
%to filter one of two excitatory inputs from reaching the soma. The
%excitatory and inhibitory inputs are square waves, of intensity
%"stimLevel" and "inhibitionIntensity." To quantify the ability to filter,
%two different metrics are used: 1) Ron's method, how does the maximum
%voltage with inhibition compare to only exciting spine1 or spine2, and
%2) Aaron's method, how does the magnitude of the inhibitory effect (still 
%exciting both spines) compare to the magnitude that exciting spine1 or 
%spine2 individually contribute. Both are passed through a raised cosine
%filter, although another bell-shaped filter (eg Gaussian) or even root
%mean squared could work.

%A cosine filter is used in order to punish "under" or "over" inhibitory
%effects. Ie, say spine1 contributes much more than spine2, then if you 
%were to inhibit spine1 then the soma voltage goes down by a lot. 
%Therefore, the max soma voltage with exc spine1/2 inh1 would be less than 
%exc spine1, and without passing it through a filter the math would suggest
%"overeffective" filtering, ie you would get a value via Ron's method less 
%than 1 (all should be >1) or via Aaron's method of sometimes >10 or more 
%(all should be between 1 and ~3).

%%

clear;
% close all;

repeats = 10;


%% Load dendritic arbor parameters
[miscParams dendParams condParams ] = loadParameters;

%% Design simulation parameters
stimLevel = 6;
inhibitionIntensity = 20; %was 5

pathCompartments = cell(repeats,2); %1 is path1, 2 is path2
path1Length = zeros(1,repeats); %holds the length of each repeat's (arbor's) path length, since it's always from 0 to 1
path1SomaDiff = cell(repeats,3); %effect of inhibition at each of the compartments along the path
path2Length = zeros(1,repeats);
path2SomaDiff = cell(repeats,3);
branchCompartment = zeros(1,repeats); %the compartment at which the two paths meet


somaMax = zeros(repeats,5); %path1, path2, path1/2, path1/2 inh1, path1/2 inh2
somaAUC = zeros(repeats,5); %path1, path2, path1/2, path1/2 inh1, path1/2 inh2
filtIndexRon = zeros(repeats,6);
filtIndexAaron = cell(repeats,2);

spineList = zeros(repeats,2);
for r = 1:repeats
    
    r
    
    %% Find spines to excite
    potentialSpinePair = [];
    
    while isempty(potentialSpinePair)
        r
        %% Build arbor
        [connectome compartmentIDs conductanceMat distance] = buildDendriticArbor(dendParams);
        g = graph(connectome);
        graphMaster(r).arbor = g;
        totalCompartments = length(compartmentIDs);
        shaftCmpt = compartmentIDs(1,find(compartmentIDs(2,:)~=4));
        spineCmpt = compartmentIDs(1,find(compartmentIDs(2,:)==4));
        
        
        %Find potential spine pairs
        for s1 = 1:length(spineCmpt)
            spineA = spineCmpt(s1);
            
            for s2 = 1:length(spineCmpt)
                spineB = spineCmpt(s2);
                
                pathA = shortestpath(g,spineA,1);pathA = pathA(2:end);
                pathB = shortestpath(g,spineB,1);pathB = pathB(2:end);
                
                if min(ismember(pathA,pathB))~=1 && min(ismember(pathB,pathA))~=1 && max(intersect(pathA,pathB))~=1
                    potentialSpinePair = [potentialSpinePair; spineA spineB];
                end
            end
        end
    end
    tempCmpt = potentialSpinePair(randi(length(potentialSpinePair)),:);
    spine1 = tempCmpt(1);
    spine2 = tempCmpt(2);
    path1 = shortestpath(g,spine1,1);path1 = path1(2:end);
    path2 = shortestpath(g,spine2,1);path2 = path2(2:end);
    spineList(r,:) = [spine1 spine2];
    shaft1 = find(connectome(spine1,:)==1);
    shaft2 = find(connectome(spine2,:)==1);
    branchCompartment(1,r) = max(intersect(path1,path2));
    pathCompartments{r,1} = path1;
    pathCompartments{r,2} = path2;
    path1Length(r) = length(path1);
    path2Length(r) = length(path2);
    
   
    %% Design stimulus
    spine1Stimulus = zeros(1,miscParams.time);
    spine1Stimulus(1,150:200) = stimLevel;
    spine2Stimulus = zeros(1,miscParams.time);
    spine2Stimulus(1,150:200) = stimLevel;
    
    %% Define pure excitatory dynamics
    iExcite = zeros(totalCompartments,miscParams.time);
    iExcite(spine1,:) = spine1Stimulus;
    iInhibit = zeros(totalCompartments,miscParams.time);
    input.excitation = iExcite;
    input.inhibition = iInhibit;
    voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
    v1 = voltages(1,:);
    somaMax(r,1) = max(v1)+60;
    somaAUC(r,1) = sum(v1+60);
    
    iExcite = zeros(totalCompartments,miscParams.time);
    iExcite(spine2,:) = spine2Stimulus;
    iInhibit = zeros(totalCompartments,miscParams.time);
    input.excitation = iExcite;
    input.inhibition = iInhibit;
    voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
    v2 = voltages(1,:);
    somaMax(r,2) = max(v2)+60;
    somaAUC(r,2) = sum(v2+60);
    
    iExcite = zeros(totalCompartments,miscParams.time);
    iExcite(spine1,:) = spine1Stimulus;
    iExcite(spine2,:) = spine2Stimulus;
    iInhibit = zeros(totalCompartments,miscParams.time);
    input.excitation = iExcite;
    input.inhibition = iInhibit;
    voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
    v12 = voltages(1,:);
    somaMax(r,3) = max(v12)+60;
    somaAUC(r,3) = sum(v12+60);
    
    % figure;hold on;
    % plot(v1); plot(v2); plot(v12);
    % legend({'V1','V2','V12'});
    
    %% Sweep through shafts to inhibit
    
    %Last dimension "3" because excite spine1 (1), spine2 (2), or both (3)
    somaTraces = zeros(length(shaftCmpt),miscParams.time,3);tempFiltIndex = zeros(2,length(shaftCmpt));
    tempFiltAaronCos = zeros(2,length(shaftCmpt));
    
    %instead of sweeping through all shafts, just sweep through the shafts
    %of path1 and path2
%     for s = 1:length(shaftCmpt) %sweep through all shafts
    u = union(path1,path2);
    for sh = 1:length(u) %sweep through shafts of path1 and path2 (u for union)
        s = u(sh);
        shaft = shaftCmpt(s);
        
        %Excite spine 1
        iExcite = zeros(totalCompartments,miscParams.time);
        iExcite(spine1,:) = spine1Stimulus;
        iInhibit = zeros(totalCompartments,miscParams.time);
        iInhibit(shaft,150:200) = inhibitionIntensity;
        input.excitation = iExcite;
        input.inhibition = iInhibit;
        voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
        
        somaTraces(s,:,1) = voltages(1,:);
        
        %Excite spine 2
        iExcite = zeros(totalCompartments,miscParams.time);
        iExcite(spine2,:) = spine2Stimulus;
        iInhibit = zeros(totalCompartments,miscParams.time);
        iInhibit(shaft,150:200) = inhibitionIntensity;
        input.excitation = iExcite;
        input.inhibition = iInhibit;
        voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
        
        somaTraces(s,:,2) = voltages(1,:);
        
        %Excite spines 1 and 2
        iExcite = zeros(totalCompartments,miscParams.time);
        iExcite(spine1,:) = spine1Stimulus;
        iExcite(spine2,:) = spine2Stimulus;
        iInhibit = zeros(totalCompartments,miscParams.time);
        iInhibit(shaft,150:200) = inhibitionIntensity;
        input.excitation = iExcite;
        input.inhibition = iInhibit;
        voltages = runSimulation(conductanceMat,compartmentIDs,input,condParams);
        
        somaTraces(s,:,3) = voltages(1,:);

        tempFiltAaronCos(1,s) = raisedCosine(max(v12)-max(voltages(1,:)),max(v12)-max(v1),max(v12)-max(v1));
        tempFiltAaronCos(2,s) = raisedCosine(max(v12)-max(voltages(1,:)),max(v12)-max(v2),max(v12)-max(v2));
    end
    
    %% Compile and store data
    
    %Spine 1
    inhDifference = max(v1) - max(somaTraces,[],2);
    path1SomaDiff{r,1} = inhDifference(path1,1);
    path2SomaDiff{r,1} = inhDifference(path2,1);
    
    %Spine 2
    inhDifference = max(v2) - max(somaTraces,[],2);
    path1SomaDiff{r,2} = inhDifference(path1,2);
    path2SomaDiff{r,2} = inhDifference(path2,2);
    
    %Spines 1 and 2
    inhDifference = max(v12) - max(somaTraces,[],2);
    path1SomaDiff{r,3} = inhDifference(path1,3);
    path2SomaDiff{r,3} = inhDifference(path2,3);
    
    %Data about filtering  
    filtIndexAaron{r,1} = tempFiltAaronCos(:,path1);
    filtIndexAaron{r,2} = tempFiltAaronCos(:,path2);
    
    somaMax(r,4) = max(somaTraces(shaft1,:,3))+60;
    somaAUC(r,4) = sum(somaTraces(shaft1,:,3)+60);
    somaMax(r,5) = max(somaTraces(shaft2,:,3))+60;
    somaAUC(r,5) = sum(somaTraces(shaft2,:,3)+60);
    
    filtIndexRon(r,1) = raisedCosine(somaMax(r,3) / somaMax(r,2),1,somaMax(r,3)/somaMax(r,2)-1);
    filtIndexRon(r,2) = raisedCosine(somaMax(r,3) / somaMax(r,1),1,somaMax(r,3)/somaMax(r,1)-1);
    filtIndexRon(r,3) = raisedCosine(somaMax(r,4) / somaMax(r,2),1,somaMax(r,3)/somaMax(r,2)-1);
    filtIndexRon(r,4) = raisedCosine(somaMax(r,5) / somaMax(r,1),1,somaMax(r,3)/somaMax(r,1)-1);
    filtIndexRon(r,5) = raisedCosine(somaMax(r,4) / somaMax(r,1),1,somaMax(r,3)/somaMax(r,1)-1);
    filtIndexRon(r,6) = raisedCosine(somaMax(r,5) / somaMax(r,2),1,somaMax(r,3)/somaMax(r,2)-1);
end

%% Analyze and plot data

%Effect of inhibition with excitation on IPSILATERAL spine
interpResolution = 101;
xInterp = linspace(0,1,interpResolution);
path1YInterp = zeros(repeats,interpResolution);
path2YInterp = zeros(repeats,interpResolution);
for r = 1:repeats
    x = linspace(0,1,path1Length(r));
    newY1 = interp1(x,path1SomaDiff{r,1},xInterp);
    path1YInterp(r,:) = newY1;
    
    x = linspace(0,1,path2Length(r));
    newY2 = interp1(x,path2SomaDiff{r,2},xInterp);
    path2YInterp(r,:) = newY2;
end
meanPath1Y = mean(path1YInterp);
meanPath2Y = mean(path2YInterp);
figure;hold on;
plot(meanPath1Y);
plot(meanPath2Y);
legend('Path 1 - excitation on spine 1','Path 2 - excitation on spine 2');
xticks([1 interpResolution]);
xticklabels({'Distal shaft', 'Soma'});
ylabel('Effect of inhibition (mV)');
title('Effect of inhibition location with ipsilateral spine excitation');

%Effect of inhibition with excitation on CONTRALATERAL spine
interpResolution = 101;
xInterp = linspace(0,1,interpResolution);
path1YInterp = zeros(repeats,interpResolution);
path2YInterp = zeros(repeats,interpResolution);
for r = 1:repeats
    x = linspace(0,1,path1Length(r));
    newY1 = interp1(x,path1SomaDiff{r,2},xInterp);
    path1YInterp(r,:) = newY1;
    
    x = linspace(0,1,path2Length(r));
    newY2 = interp1(x,path2SomaDiff{r,1},xInterp);
    path2YInterp(r,:) = newY2;
end
meanPath1Y = mean(path1YInterp);
meanPath2Y = mean(path2YInterp);
figure;hold on;
plot(meanPath1Y);
plot(meanPath2Y);
legend('Path 1 - excitation on spine 2','Path 2 - excitation on spine 1');
xticks([1 interpResolution]);
xticklabels({'Distal shaft', 'Soma'});
ylabel('Effect of inhibition (mV)');
title('Effect of inhibition location with contralateral spine excitation');

%Effect of inhibition with excitation on IPSI/CONTRALATERAL spine
interpResolution = 101;
xInterp = linspace(0,1,interpResolution);
path1YInterp = zeros(repeats,interpResolution*2-1);
path2YInterp = zeros(repeats,interpResolution*2-1);
for r = 1:repeats
    if branchCompartment(r)~=1
        prePath = pathCompartments{r,1}(1:find(pathCompartments{r,1}==branchCompartment(r)));
        postPath = pathCompartments{r,1}(find(pathCompartments{r,1}==branchCompartment(r)):end);
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,path1SomaDiff{r,3}(1:length(prePath)),xInterp);
        x = linspace(0,1,length(postPath));
        postInterp = interp1(x,path1SomaDiff{r,3}(end-length(postPath)+1:end),xInterp);
        path1YInterp(r,:) = [preInterp postInterp(2:end)];
        
        prePath = pathCompartments{r,2}(1:find(pathCompartments{r,2}==branchCompartment(r)));
        postPath = pathCompartments{r,2}(find(pathCompartments{r,2}==branchCompartment(r)):end);
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,path2SomaDiff{r,3}(1:length(prePath)),xInterp);
        x = linspace(0,1,length(postPath));
        postInterp = interp1(x,path2SomaDiff{r,3}(end-length(postPath)+1:end),xInterp);
        path2YInterp(r,:) = [preInterp postInterp(2:end)];
    else
        prePath = pathCompartments{r,1}(1:find(pathCompartments{r,1}==branchCompartment(r)));
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,path1SomaDiff{r,3}(1:length(prePath)),xInterp);
        postInterp = ones(1,interpResolution)*path1SomaDiff{r,3}(end);
        path1YInterp(r,:) = [preInterp postInterp(2:end)];
        
        prePath = pathCompartments{r,2}(1:find(pathCompartments{r,2}==branchCompartment(r)));
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,path2SomaDiff{r,3}(1:length(prePath)),xInterp);
        postInterp = ones(1,interpResolution)*path2SomaDiff{r,3}(end);
        path2YInterp(r,:) = [preInterp postInterp(2:end)];
    end
end
meanPath1Y = mean(path1YInterp);
meanPath2Y = mean(path2YInterp);
figure;hold on;
plot(meanPath1Y);
plot(meanPath2Y);
line([101 101],[0 max([meanPath2Y meanPath1Y])]);
legend('Path 1 - excitation on spine 1/2','Path 2 - excitation on spine 1/2');
xticks([1 interpResolution interpResolution*2-1]);
xticklabels({'Distal shaft', 'Branch point','Soma'});
ylabel('Effect of inhibition (mV)');
title('Effect of inhibition location with ipsi- and contralateral spine excitation');


meanFiltIndexRonCos = mean(filtIndexRon);
figure;bar([meanFiltIndexRonCos(2) meanFiltIndexRonCos(5) meanFiltIndexRonCos(4); meanFiltIndexRonCos(1) meanFiltIndexRonCos(3) meanFiltIndexRonCos(6)]);
xticklabels({'Exc 1/2 vs 1, ±inh','Exc 1/2 vs 2, ±inh'});
legend({'No inhibition','Inhibit path 1','Inhibit path 2'});
title('Ron filtration index');

%Filtering index, looks like v1
interpResolution = 101;
xInterp = linspace(0,1,interpResolution);
path1YInterp = zeros(repeats,interpResolution*2-1);
path2YInterp = zeros(repeats,interpResolution*2-1);
for r = 1:repeats
    if branchCompartment(r)~=1
        prePath = pathCompartments{r,1}(1:find(pathCompartments{r,1}==branchCompartment(r)));
        postPath = pathCompartments{r,1}(find(pathCompartments{r,1}==branchCompartment(r)):end);
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,filtIndexAaron{r,1}(1,1:length(prePath)),xInterp);
        x = linspace(0,1,length(postPath));
        postInterp = interp1(x,filtIndexAaron{r,1}(1,end-length(postPath)+1:end),xInterp);
        path1YInterp(r,:) = [preInterp postInterp(2:end)];
        
        prePath = pathCompartments{r,2}(1:find(pathCompartments{r,2}==branchCompartment(r)));
        postPath = pathCompartments{r,2}(find(pathCompartments{r,2}==branchCompartment(r)):end);
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,filtIndexAaron{r,2}(1,1:length(prePath)),xInterp);
        x = linspace(0,1,length(postPath));
        postInterp = interp1(x,filtIndexAaron{r,2}(1,end-length(postPath)+1:end),xInterp);
        path2YInterp(r,:) = [preInterp postInterp(2:end)];
    else
        prePath = pathCompartments{r,1}(1:find(pathCompartments{r,1}==branchCompartment(r)));
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,filtIndexAaron{r,1}(1,1:length(prePath)),xInterp);
        postInterp = ones(1,interpResolution)*filtIndexAaron{r,1}(1,end);
        path1YInterp(r,:) = [preInterp postInterp(2:end)];
        
        prePath = pathCompartments{r,2}(1:find(pathCompartments{r,2}==branchCompartment(r)));
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,filtIndexAaron{r,2}(1,1:length(prePath)),xInterp);
        postInterp = ones(1,interpResolution)*filtIndexAaron{r,2}(1,end);
        path2YInterp(r,:) = [preInterp postInterp(2:end)];
    end
end
meanPath1Y = mean(path1YInterp);
meanPath2Y = mean(path2YInterp);
figure;hold on;
plot(meanPath1Y);
plot(meanPath2Y);
line([101 101],[0 max([meanPath2Y meanPath1Y])]);
legend('Path 1','Path 2');
xticks([1 interpResolution interpResolution*2-1]);
xticklabels({'Distal shaft', 'Branch point','Soma'});
ylabel('Ability to filter out spine 2');
title('Aaron filtration index - filter spine 2');

%Filtering index: looks like v2
interpResolution = 101;
xInterp = linspace(0,1,interpResolution);
path1YInterp = zeros(repeats,interpResolution*2-1);
path2YInterp = zeros(repeats,interpResolution*2-1);
for r = 1:repeats
    if branchCompartment(r)~=1
        prePath = pathCompartments{r,1}(1:find(pathCompartments{r,1}==branchCompartment(r)));
        postPath = pathCompartments{r,1}(find(pathCompartments{r,1}==branchCompartment(r)):end);
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,filtIndexAaron{r,1}(2,1:length(prePath)),xInterp);
        x = linspace(0,1,length(postPath));
        postInterp = interp1(x,filtIndexAaron{r,1}(2,end-length(postPath)+1:end),xInterp);
        path1YInterp(r,:) = [preInterp postInterp(2:end)];
        
        prePath = pathCompartments{r,2}(1:find(pathCompartments{r,2}==branchCompartment(r)));
        postPath = pathCompartments{r,2}(find(pathCompartments{r,2}==branchCompartment(r)):end);
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,filtIndexAaron{r,2}(2,1:length(prePath)),xInterp);
        x = linspace(0,1,length(postPath));
        postInterp = interp1(x,filtIndexAaron{r,2}(2,end-length(postPath)+1:end),xInterp);
        path2YInterp(r,:) = [preInterp postInterp(2:end)];
    else
        prePath = pathCompartments{r,1}(1:find(pathCompartments{r,1}==branchCompartment(r)));
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,filtIndexAaron{r,1}(2,1:length(prePath)),xInterp);
        postInterp = ones(1,interpResolution)*filtIndexAaron{r,1}(2,end);
        path1YInterp(r,:) = [preInterp postInterp(2:end)];
        
        prePath = pathCompartments{r,2}(1:find(pathCompartments{r,2}==branchCompartment(r)));
        x = linspace(0,1,length(prePath));
        preInterp = interp1(x,filtIndexAaron{r,2}(2,1:length(prePath)),xInterp);
        postInterp = ones(1,interpResolution)*filtIndexAaron{r,2}(2,end);
        path2YInterp(r,:) = [preInterp postInterp(2:end)];
    end
end
meanPath1Y = mean(path1YInterp);
meanPath2Y = mean(path2YInterp);
figure;hold on;
plot(meanPath1Y);
plot(meanPath2Y);
line([101 101],[0 max([meanPath2Y meanPath1Y])]);
legend('Path 1','Path 2');
xticks([1 interpResolution interpResolution*2-1]);
xticklabels({'Distal shaft', 'Branch point','Soma'});
ylabel('Ability to filter out spine 1');
title('Aaron filtration index - filter spine 1');


function output = raisedCosine(input,cosU,inS)
    cosS = inS;
    if input<(cosU+cosS) && input>(cosU-cosS)
        output = 1/(2) * (1+ cos((input-cosU)/cosS*pi)); %1/2 NOT 1/2s so that max is 1
    else
        output = 0;
    end
end