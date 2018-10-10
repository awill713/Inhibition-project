function [ voltage ] = runSimulation( dendConduct, compIDs, extInput, par)
%RUNSIMULATION Summary of this function goes here
%   Faster implementation, clean variables naming/code later

conductances = dendConduct;
compIdentity = compIDs;
excite = extInput.excitation;
inhibit = extInput.inhibition;
params = par;
p = params; %fix later
totalCompartments = size(compIdentity,2);

edges = [1 zeros(1, max(compIdentity(2,:)+1)-2) totalCompartments];
for inds = 2 : length(edges)-1
    
    edges(inds) = find(compIdentity(2,:) == inds-1, 1, 'last');
    
end

group_masks = zeros(totalCompartments,1, max(compIdentity(2,:)+1));
group_masks(1,1,1)= 1;

for inds = 2: max(compIdentity(2,:)+1)
    
    group_masks(edges(inds-1)+1:edges(inds),1,inds) = ones(diff(edges(inds-1:inds)),1);
end


voltage = zeros(totalCompartments,params.timePoints);
voltage(:,1) = par.vRest;
u = zeros(totalCompartments,params.timePoints);


%alternatively could do a mask approach based on edge
%Now doing this

for t = 2:params.timePoints
  
  pvt = (p.shaftK * ( voltage(:,t-1)  - p.vRest).*group_masks(:,:,1) .* (voltage(:,t-1)-p.vThresh).* group_masks(:,:,1) - u(:,t-1).*group_masks(:,:,1))+...
       p.shaftK * ( p.vRest - voltage(:,t-1)).*group_masks(:,:,2) +...
       p.shaftK * ( p.vRest - voltage(:,t-1)).*group_masks(:,:,3) +...
       p.shaftK * ( p.vRest - voltage(:,t-1)).*group_masks(:,:,4) +...
      (p.spineK * ( voltage(:,t-1)  - p.vRest).*group_masks(:,:,5) .* (voltage(:,t-1)-p.vThresh).* group_masks(:,:,5) - u(:,t-1).*group_masks(:,:,5));
      
  
  
  input_term = excite(:,t-1).*(p.vNa - voltage(:,t-1)) + inhibit(:,t-1).*(p.vCl - voltage(:,t-1));
  
  neighbor_term = sum(conductances.*(meshgrid(voltage(:,t-1)) - meshgrid(voltage(:,t-1))'),2);
  
  dv = ((pvt + input_term + neighbor_term).*group_masks(:,:,1)*p.dt/p.C) +...soma has to be done differently
      ((pvt + input_term + neighbor_term) .*sum(group_masks(:,:,2:5),3)*p.dt/p.Celse);
  
  voltage(:,t) =  voltage(:,t-1) + dv;
  
  du = p.a*(p.b*(voltage(:,t-1) - p.vRest) - u(:,t-1)) * p.dt;
  u(:,t) = u(:,t-1) + du;
  
  u((voltage(:,t)>p.vPeak),t) = u((voltage(:,t)>p.vPeak),t) + p.d;
  voltage((voltage(:,t)>p.vPeak),t) = p.c;
  

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
