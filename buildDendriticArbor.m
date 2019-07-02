function [ connect cIDs conduct dist] = buildDendriticArbor( input )
%BUILDDENDRITICARBOR Summary of this function goes here
%   Detailed explanation goes here

parameters = input;

connect = 0;
cIDs = [1;0];
conduct = 0;

%% Basal dendrites
level = 0;
compType = 1;

for bBranch = 1:parameters.basalBranches
    
    level = level + 1;
    tempParents = [];
    
    [connect cIDs tempParents conduct] = addCompartment(connect, cIDs, compType,1,tempParents,conduct,parameters);
    
    parentCompartments = tempParents;
    
    while level < parameters.basalDepth && length(parentCompartments)>0
        
        level = level + 1;
        tempParents = [];
        
        for p = 1:length(parentCompartments)
            parentID = parentCompartments(p);
            
            daughters = randsample(0:1:parameters.basalDaughters,1,true,parameters.basalProb);
            
            for d = 1:daughters
                [connect cIDs tempParents conduct] = addCompartment(connect,cIDs,compType,parentID,tempParents,conduct,parameters);
            end
        end
        
        parentCompartments = tempParents;
    end
    
    level = 0;
    
end


%% Apical shaft
root = 1;
compType = 2;
for a = 1:parameters.apicalShaftCompartments
    level = 0;
    tempParents = [];
    
    [connect cIDs tempParents conduct] = addCompartment(connect,cIDs,compType,root,tempParents,conduct,parameters);

    root = tempParents; %tempParents happens to be just this compartment
    parentCompartments = [tempParents];

    while level < parameters.apicalDepth && length(parentCompartments)>0
        
        level = level + 1;
        tempParents = [];
        
        for p = 1:length(parentCompartments)
            parentID = parentCompartments(p);
            
            daughters = randsample(0:1:parameters.apicalDaughters,1,true,parameters.apicalProb);
%                 [num2str(parentID) ' has ' num2str(daughters) ' daughters']
            for d = 1:daughters
                [connect cIDs tempParents conduct] = addCompartment(connect,cIDs,compType,parentID,tempParents,conduct,parameters);
            end
        end
        
        parentCompartments = tempParents;
    end
    
    level = 0;
    
end

%% Apical tuft
level = 0;
parentCompartments = root;
compType = 3;
% parentCompartments
% parameters.tuftDepth
% level
while level < parameters.tuftDepth && length(parentCompartments)>0
    
    level = level + 1;
    tempParents = [];
    
    for p = 1:length(parentCompartments)
        parentID = parentCompartments(p);
        
        daughters = randsample(0:1:parameters.tuftDaughters,1,true,parameters.tuftProb);
%         [num2str(parentID) ' has ' num2str(daughters) ' daughters']
        for d = 1:daughters
            [connect cIDs tempParents conduct] = addCompartment(connect,cIDs,compType,parentID,tempParents,conduct,parameters);
        end
    end
    
    parentCompartments = tempParents;
%     display('Did apical tuft');
end

%% Add dendritic spines
compType = 4;
for c = 2:size(cIDs,2)
    for s = 1:parameters.spinesPerCompartment
        [connect cIDs tempParents conduct] = addCompartment(connect,cIDs,compType,c,tempParents,conduct,parameters);
    end
end

g = graph(connect);
dist = distances(g);

%% Visualize connectome structure
% g = graph(connect);
% figure; h = plot(g);
% highlight(h,cIDs(1,find(cIDs(2,:)==1)),'NodeColor','g');
% highlight(h,cIDs(1,find(cIDs(2,:)==2)),'NodeColor','k');
% highlight(h,cIDs(1,find(cIDs(2,:)==3)),'NodeColor','r');
% title('Green is basal dendrites, black is apical, red is apical tuft');

end


function [outC outID outTP outConduct] = addCompartment(c, id, type, parent, tp, con, para)
    outC = c;
    outID = id;
    outTP = tp;
    outConduct = con;
    
    currentID = length(outID(1,:)) + 1;
    outID = [outID [currentID;type]];
    outC(parent,currentID) = 1;
    outC(currentID,parent) = 1;
    outTP = [outTP currentID];
    
    if type == 4
        g = para.spineConduct;
    else
        g = para.shaftConduct;
    end
    
    outConduct(parent,currentID) = g;
    outConduct(currentID,parent) = g;
end


