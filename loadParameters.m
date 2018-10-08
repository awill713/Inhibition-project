function [ misc dend cond ] = loadParameters( input_args )
%LOADPARAMETERS Summary of this function goes here
%   Detailed explanation goes here

%% Miscellaneous parameters
misc.totalTime=50;%length of sim in ms
misc.dt = 0.1; %ms
misc.time = round(misc.totalTime/misc.dt);

%% Dendritic tree initialization parameters

dend.basalBranches = 3; %roots branching directly from the soma
dend.basalDepth = 3; %maximum number of compartments away from soma
dend.basalDaughters = 2; %how many daughter branches at each branch point
dend.basalProb = [0.1 0.4 0.5]; %probability of number of branches at each branch point (should have length basalDaughters+1)

dend.apicalShaftCompartments = 5; %length of apical dendrite in compartments
dend.apicalDepth = 2; %maximum number of compartments away from main apical shaft
dend.apicalBranches = 2; %roots branching directly from each apical shaft compartment
dend.apicalDaughters = 2; %how many daughter branches at each branch point
dend.apicalProb = [0.7 0.25 0.05]; %probability of number of branches at each branch point (should have length apicalDaughters+1)

dend.tuftDepth = 5; %maximum number of compartments away from last apical shaft compartment
dend.tuftDaughters = 2; %how many daughter branches at each branch point

%changed this variable from [0.05 0.6 0.35] because the 5% chance of no
%apical shaft compartments threw error with runSimulation edges
dend.tuftProb = [0.00 0.6 0.4]; %probability of number of branches at each branch point (should have length tuftDaughters+1)

dend.spinesPerCompartment = 1; %number of spines connected to each compartment, except for the soma
dend.shaftConduct = 150; %conductance between shaft (and soma) compartments
dend.spineConduct = 25; %conductance between spine and shaft compartments

%% Intercompartment dynamics parameters
cond.C = 100;
cond.Celse = 50;
cond.vRest = -60;
cond.vThresh = -50;
cond.a = 0.01;
cond.b = -2;
cond.c = -55;
cond.d = 100;
cond.vPeak = 50;
cond.vNa = 0;
cond.vCl = -70;
cond.shaftK = 0.7;
cond.spineK = 3;
% conductParams.shaftConduct = 150;
% conductParams.spineConduct = 25;
cond.dt = misc.dt;
cond.timePoints = misc.time;

end

