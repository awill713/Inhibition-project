
connectome = 0;
compartmentIDs = [1;0];

%% Basal dendrites
basalBranches = 3;
basalDepth = 3;
basalDaughters = 2;
basalProb = [0.1 0.4 0.5];
% basalProb = [1 0 0];

level = 0;

for bBranch = 1:basalBranches
    
    level = level + 1;
    
    currentID = length(compartmentIDs(1,:))+1;
    compartmentIDs = [compartmentIDs [currentID;1]];
    connectome(1,currentID) = 1;
    connectome(currentID,1) = 1;
    parentCompartments = [currentID];
    
    while level < basalDepth && length(parentCompartments)>0
        
        level = level + 1;
        tempParents = [];
        
        for p = 1:length(parentCompartments)
            parentID = parentCompartments(p);
            
            daughters = randsample(0:1:basalDaughters,1,true,basalProb);
            
            for d = 1:daughters
                currentID = length(compartmentIDs(1,:)) + 1;
                compartmentIDs = [compartmentIDs [currentID;1]];
                connectome(parentID,currentID) = 1;
                connectome(currentID,parentID) = 1;
                tempParents = [tempParents currentID];
            end
        end
        
        parentCompartments = tempParents;
    end
    
    level = 0;
    
end


%% Apical shaft

apicalShaftCompartments = 5;
apicalDepth = 2;
apicalBranches = 2; %per shaft compartment
apicalDaughters = 2;
apicalProb = [0.7 0.25 0.05];
% apicalProb = [1 0 0];

root = 1;
for a = 1:apicalShaftCompartments
    currentID = length(compartmentIDs(1,:))+1;
    compartmentIDs = [compartmentIDs [currentID;2]];
    connectome(root,currentID) = 1;
    connectome(currentID,root) = 1;
    root = currentID;
    parentCompartments = [currentID];
    
    level = 0;
    
%     for aBranch = 1:apicalBranches
%         aBranch
%         
%         level = level + 1;
        
%         currentID = length(compartmentIDs(1,:))+1;
%         compartmentIDs = [compartmentIDs [currentID;2]];
%         connectome(root,currentID) = 1;
%         connectome(currentID,root) = 1;
%         parentCompartments = [currentID];
        
        while level < apicalDepth && length(parentCompartments)>0
            
            level = level + 1;
            tempParents = [];
            
            for p = 1:length(parentCompartments)
                parentID = parentCompartments(p);
                
                daughters = randsample(0:1:apicalDaughters,1,true,apicalProb);
                ['Level ' num2str(level) ' c ' num2str(parentID) ' has ' num2str(daughters) ' daughters']
                for d = 1:daughters
                    currentID = length(compartmentIDs(1,:)) + 1;
                    compartmentIDs = [compartmentIDs [currentID;2]];
                    connectome(parentID,currentID) = 1;
                    connectome(currentID,parentID) = 1;
                    tempParents = [tempParents currentID];
                end
            end
            
            parentCompartments = tempParents;
            level
        end
        
        level = 0
        
%     end
end

%% Apical tuft
tuftDepth = 5;
tuftDaughters = 2;
tuftProb = [0.05 0.6 0.35];
% tuftProb = [1 0 0];

level = 0;
parentCompartments = root;

while level < tuftDepth && length(parentCompartments)>0
    
    level = level + 1;
    tempParents = [];
    
    for p = 1:length(parentCompartments)
        parentID = parentCompartments(p);
        
        daughters = randsample(0:1:tuftDaughters,1,true,tuftProb);
%         [num2str(parentID) ' has ' num2str(daughters) ' daughters']
        for d = 1:daughters
            currentID = length(compartmentIDs(1,:)) + 1;
            compartmentIDs = [compartmentIDs [currentID;3]];
            connectome(parentID,currentID) = 1;
            connectome(currentID,parentID) = 1;
            tempParents = [tempParents currentID];
        end
    end
    
    parentCompartments = tempParents;
end

%% Visualize connectome structure
g = graph(connectome);
figure; h = plot(g);
highlight(h,compartmentIDs(1,find(compartmentIDs(2,:)==1)),'NodeColor','g');
highlight(h,compartmentIDs(1,find(compartmentIDs(2,:)==2)),'NodeColor','k');
highlight(h,compartmentIDs(1,find(compartmentIDs(2,:)==3)),'NodeColor','r');

%% Add dendritic spines
spinesPerShaft = 1;

for c = 2:size(compartmentIDs,2)
    for s = 1:spinesPerShaft
        currentID = length(compartmentIDs(1,:)) + 1;
        compartmentIDs = [compartmentIDs [currentID;4]];
        connectome(c,currentID) = 1;
        connectome(currentID,c) = 1;
    end
end