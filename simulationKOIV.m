%% SIMULATION SCRIPT FOR CALLING KOIVUMAKI MODEL
% Outputs with the time, state variables and currents matrices for a simulation.
% Natalia Trigueros
function [time,X,IsJs] = simulationKOIV(scaling_factors, BeatsSaved, pAF_state)

% This function simulates the action potential model
% 60-beat simulation with a basic cycle length of 1000 ms
beats=60; BCL=1; 

% Create matrices for storing data retrieved from model
time = cell(1,BeatsSaved);
X = cell(1,BeatsSaved);
IsJs= cell(1,BeatsSaved); % Currents storage
StateVars_Court= cell(1,BeatsSaved); % State variables storage
Ti_Court=cell(1,BeatsSaved);  % Time storage

% Import initial conditions for Koivumaki 2011 model from .txt
y0_filename = 'y0_ss_BCL_1000_koivumaki2011_TableS4.dat.txt';
x0 = importdata(y0_filename);

% Set options
options = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'MaxStep', 0.9e-3); 
ACh = "ACh005"; % Go to model to see options for ACh and mutation
mutation = 'NO';
pos=1;

for n=1:beats % n = 1:60
    if n <= beats-BeatsSaved
        [time_null, X_null]= ode15s(@Koivumaki_model,[0 1], x0, options, ACh, mutation, scaling_factors, pAF_state);
        x0=X_null(end,:);
        
    elseif n > beats-BeatsSaved % Enter this iteration beats 58, 59 and 60
        [time{1,pos}, X{1,pos}] = ode15s(@Koivumaki_model, [0 1], x0, options, ACh, mutation, scaling_factors, pAF_state);
        x0=X{1,pos}(size(X{1,pos},1),:);
        pos=pos+1;
    end
end

for b=1:BeatsSaved % b = 1:3
    for i=1:size(X{b},1)
        [StateVars_Court{1,b}(i,:),IsJs{1,b}(i,:),Ti_Court{1,b}(i,:)] = Koivumaki_model(time{b}(i),X{b}(i,:), ACh, mutation, scaling_factors, pAF_state);
        % Retrieve from model state vars, currents and time
    end
end

end