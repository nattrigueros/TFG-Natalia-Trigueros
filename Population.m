%% MAIN FOR POPULATION CREATION: CRN (Control/pAF), KOIVUMAKI(Control/pAF)
% Natalia Trigueros

trials = 1000; % number of trials (rows); and parameters (columns)
%% STEP 1: Choose model and remodelling option through command window.

m_n = input('Choose model from options: a / b \n\na.CRN \nb.KOIVUMAKI \n\nYour selection is:','s');

switch m_n
    case {'a','a.','A','A.','CRN','a.CRN'}
        model_name='CRN'
    case {'b','b.','B','B.','KOIVUMAKI','b.KOIVUMAKI'}
        model_name='KOIVUMAKI'
    otherwise 
        disp(['Run again file and choose a correct option a/b'])
end

pAF_state='pAF_no'

%% STEP 2: Parameter factors generation.
% Create scaling factors matrix that will follow a normal distribution.
% Scaling factors matrix will have as many columns as parameter in output
% (currents)of model.

switch model_name
    
    case 'CRN'
        parameter_names = {'i_Na', 'i_K1', 'i_to', 'i_Kur', 'i_Kr', 'i_Ks', 'i_Ca_L', 'i_NaK', 'i_NaCa', 'i_CaP',...
        'i_B_Na','i_B_Ca', 'i_B_K',...
        'i_rel', 'i_tr', 'i_up_leak', 'i_up', 'i_KAch'};
   
    case 'KOIVUMAKI'
        parameter_names={'Iion', 'Istim', 'INa', 'ICaL', 'It', 'Isus', 'IKr', 'IKs', 'IK1', 'If', 'INab', 'ICab', 'ICaP', 'INaK', 'INaCa','i_KACh'};
    
    otherwise 
        disp(['Run again file and choose a correct option a/b'])
end
%                 

n_parameters = length(parameter_names); % number of trials (rows); and parameters (columns)
mu=1; % Mean
sigma = 0.2; % Standard deviation
rng('shuffle');
p=normrnd(mu, sigma, [trials,n_parameters]); % Scaling factors matrix

%% STEP 3: Parallelization of computation.
% Throughout this file we have used paralellization. In the following lines
% we initialize workers and perform simulations in parallel (1 per worker)
% aiming a faster performance.

tic
myCluster=parcluster('local');
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool(myCluster.NumWorkers);
end
poolobj.IdleTimeout= 120;
disp(datetime)
disp(['Usando ' num2str(poolobj.NumWorkers) 'workers'])

%% STEP 4: Perform the simulations.
% Depending on model and remodelling option we will call a certain function
% and pass certain parameters.
% Simulations will be saved in .mat for further study.
switch model_name
    case 'CRN'
        parfor k=1:trials
            p_SA = p(k,:); 
            [time{k},X{k},IsJs{k}] = simulationCRN(p_SA, 3, pAF_state)
        end
        
        % Save simulation and erase unnecessary elements
        switch pAF_state
            case 'pAF_no'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationCRN_raw';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationCRN_raw.mat'], 's')
        end
    case 'KOIVUMAKI'
        parfor k=1:trials
            p_SA = p(k,:); % 17 parameters
            [time{k},X{k},IsJs{k}] = simulationKOIV(p_SA, 3, pAF_state)
        end
        
        % Save simulation and erase unnecessary elements
        switch pAF_state
            case 'pAF_no'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationKOIV_raw';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationKOIV_raw.mat'], 's')

        end
end
clearvars -except p trials
