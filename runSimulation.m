function [ voltage ] = runSimulation( dendConduct, compIDs, extInput, par)
%RUNSIMULATION Summary of this function goes here
%   Detailed explanation goes here

conductances = dendConduct;
compIdentity = compIDs;
excite = extInput.excitation;
inhibit = extInput.inhibition;
params = par;

totalCompartments = size(compIdentity,2);

voltage = zeros(totalCompartments,params.timePoints);
voltage(:,1) = par.vRest;
u = zeros(totalCompartments,params.timePoints);

for t = 2:params.timePoints
    for c = 1:totalCompartments
        [voltage(:,t) u(:,t)] = timeStep(conductances, compIdentity, voltage(:,t-1), u(:,t-1), excite(:,t-1), inhibit(:,t-1),params);
    end
end

end

function [v uu] = timeStep(conduct, cIDs, volt, you, exc, inh, p)

v = zeros(length(volt),1);
uu = zeros(length(volt),1);

for comp = 1:length(v)
    vNow = volt(comp);
    
    if cIDs(2,comp) == 4
        k = p.spineK;
        CC = p.Celse;
        previousVoltageTerm = k*(vNow-p.vRest)*(vNow-p.vThresh) - you(comp);
    elseif cIDs(2,comp) == 0
        k = p.shaftK;
        CC = p.C;
        previousVoltageTerm = k*(vNow-p.vRest)*(vNow-p.vThresh) - you(comp);
    else
        k = p.shaftK;
        CC = p.Celse;
        %         previousVoltageTerm = (vNow-p.vRest);
        previousVoltageTerm = (p.vRest-vNow);
    end
    
    inputTerm = exc(comp)*(p.vNa - vNow) + inh(comp)*(p.vCl - vNow);
    
    neighborsTerm = sum(conduct(comp,:)' .* (volt - vNow));
    
    dv = (previousVoltageTerm + inputTerm + neighborsTerm) * p.dt/CC;
    
    %     dv = (p.k*(vNow-vRest)*(vNow - vThresh) - you(comp) + sum(conduct(comp,:)'.*(volt - vNow))) * dt / C;
    v(comp) = vNow + dv;
    
    du = p.a*(p.b*(vNow - p.vRest) - you(comp)) * p.dt;
    uu(comp) = you(comp) + du;
    
    if v(comp)>p.vPeak
        v(comp) = p.c;
        uu(comp) = you(comp) + p.d;
    end
end

end
