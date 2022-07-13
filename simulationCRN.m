%% SIMULATION SCRIPT FOR CALLING CRN MODEL
% Outputs with the time, state variables and currents matrices for a 
% simulation.
% Natalia Trigueros

function [time,X,IsJs] = simulationCRN(scaling_factors, BeatsSaved, pAF_state)
% This function simulates the action potential model
% 60-beat simulation with a basic cycle length of 1000 ms

beats=60; BCL=1000; 
numStim=60; 

% Create matrices for storing data retrieved from model
time = cell(1,BeatsSaved);
X = cell(1,BeatsSaved);
IsJs= cell(1,BeatsSaved); % Currents storage
StateVars_Court= cell(1,BeatsSaved); % State variables storage
Ti_Court=cell(1,BeatsSaved);  % Time storage

% Define settings
offset=0.0; % sec
threshold=  -754;  % Threshold Courtemanche  
amplitudeStim=2*threshold; 
durationStim=2; % msec
ACh = "ACh005";
opts= odeset('RelTol',1e-6); 

% Define initial conditions for CRN model
settings = setDefaultSettings(numStim,offset,BCL,amplitudeStim,durationStim,pAF_state);
x0=getInitialVector(settings); 	% Initial conditions
pos=1;

for n=1:beats % n = 1:60
    if n <= beats-BeatsSaved
        [time_null, X_null]=ode15s('Courtemanche_model',[0 settings.BCL],x0,opts,settings,ACh,scaling_factors,pAF_state); %x0 is Y in R_courtemanche_model
        x0=X_null(end,:);
        
    elseif n > beats-BeatsSaved % Enter this iteration beats 58, 59 and 60
        [time{1,pos}, X{1,pos}]=ode15s('Courtemanche_model',[0 settings.BCL],x0,opts,settings,ACh,scaling_factors,pAF_state); %x0 is Y in R_courtemanche_model
        x0=X{1,pos}(size(X{1,pos},1),:);
        pos=pos+1;
    end
end

for b=1:BeatsSaved  % b = 1:3
    for i=1:size(X{b},1)
        [StateVars_Court{1,b}(i,:),IsJs{1,b}(i,:),Ti_Court{1,b}(i,:)] = Courtemanche_model(time{b}(i),X{b}(i,:),0,settings,ACh,scaling_factors,pAF_state);
        % Retrieve from model state vars, currents and time
    end
end

function StateInit = getInitialVector(settings)
    RMP_myo=-81.18;  
    % Initial Conditions from the paper
    StateInit = [0 1 0.999 1.37e-4 0.775 0.999 0.965 0.978 2.91e-3 1.02e-4 1.49 1.49 1.39e2 11.2 RMP_myo 3.29e-5 1.87e-2 3.04e-2 0.999 4.96e-3 0.999];
    % YNames = {'u', 'v', 'w', 'd', 'f_Ca', 'f', 'h', 'j', 'm', 'Ca_i', 'Ca_rel', 'Ca_up', 'K_i', 'Na_i', 'V', 'xr', 'xs', 'oa', 'oi', 'ua', 'ui'};

end


function settings = setDefaultSettings(numStim,offset,BCL,amplitudeStim,durationStim,pAF_state)
settings.start=1; % Initialize settings parameter
% 	if ~isfield(settings, 'dt'), settings.dt = 0.005; end		% time step	ms uncomment if using a fix step
	if ~isfield(settings, 'NumStim'), settings.NumStim= numStim; end
    if ~isfield(settings, 'storeLast'), settings.storeLast = settings.NumStim; end
	if ~isfield(settings, 'StimOffset'), settings.StimOffset = offset; end		% Offset of the first stimuls in milliseconds
    if ~isfield(settings, 'BCL'), settings.BCL = BCL; end % BCL=Base Cycle Length
    if ~isfield(settings, 'Amp_stim'), settings.Amp_stim =amplitudeStim; end % Amp stim in pA
    if ~isfield(settings, 'Dur_stim'), settings.Dur_stim = durationStim; end
    if ~isfield(settings, 'pAF'), settings.pAF = pAF_state; end %1 -47.8, 2 -31.4 mV  Resting Membrane Potential

end

end