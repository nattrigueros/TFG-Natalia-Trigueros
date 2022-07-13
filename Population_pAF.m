%% MAIN FOR POPULATION CREATION: CRN (Control/pAF), KOIVUMAKI(Control/pAF)
% Natalia Trigueros
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

region = input('Choose atrial region from options: a / b / c / d / e \n\na.PV \nb.BBLA \nc.LAA \nd.MVR \ne.LA\n\nYour selection is:','s');

switch region
    case {'a'}
        atrial_region='PV'
    case {'b'}
        atrial_region='BBLA'
    case {'c'}
        atrial_region='LAA'
    case {'d'}
        atrial_region='MVR'
    case {'e'}
        atrial_region='LA'
    otherwise 
        disp(['Run again file and choose a correct option a/b'])
end

pAF_state='pAF_yes'

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

switch model_name
    case 'CRN'
        switch atrial_region
            case 'PV'
                f=load (['ScalingFactorsCRN_PV.mat']);
                pv=f.pv;
                trials = size(pv.params,1);
                p = pv.params;
        
            case 'BBLA'
                f=load (['ScalingFactorsCRN_BBLA.mat']);
                bbla=f.bbla;
                trials = size(bbla.params,1);
                p = bbla.params;
            
            case 'LAA'
                f=load (['ScalingFactorsCRN_LAA.mat']);
                laa=f.laa;
                trials = size(laa.params,1);
                p = laa.params;
        
            case 'MVR'
                f=load (['ScalingFactorsCRN_MVR.mat']);
                mvr=f.mvr;
                trials = size(mvr.params,1);
                p = mvr.params;
        
            case 'LA'
                f=load (['ScalingFactorsCRN_LA.mat']);
                la=f.la;
                trials = size(la.params,1);
                p = la.params;

            otherwise 
                disp(['Run again file and choose a correct option a/b'])
        end

    case 'KOIVUMAKI'
        switch atrial_region
            case 'PV'
                f=load (['ScalingFactorsKOIV_PV.mat']);
                pv=f.pv;
                trials = size(pv.params,1);
                p = pv.params;
        
            case 'BBLA'
                f=load (['ScalingFactorsKOIV_BBLA.mat']);
                bbla=f.bbla;
                trials = size(bbla.params,1);
                p = bbla.params;
            
            case 'LAA'
                f=load (['ScalingFactorsKOIV_LAA.mat']);
                laa=f.laa;
                trials = size(laa.params,1);
                p = laa.params;
        
            case 'MVR'
                f=load (['ScalingFactorsKOIV_MVR.mat']);
                mvr=f.mvr;
                trials = size(mvr.params,1);
                p = mvr.params;
        
            case 'LA'
                f=load (['ScalingFactorsKOIV_LA.mat']);
                la=f.la;
                trials = size(la.params,1);
                p = la.params;

            otherwise 
                disp(['Run again file and choose a correct option a/b'])
        end
end
%                 

n_parameters = length(parameter_names); % number of trials (rows); and parameters (columns)

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
            [time{k},X{k},IsJs{k}] = simulationCRN(p(k,:), 3, pAF_state);
        end
        
        % Save simulation and erase unnecessary elements
        switch atrial_region
            case 'PV'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationCRN_pAF_PV';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationCRN_pAF_PV.mat'], 's')
            case 'BBLA'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationCRN_pAF_BBLA';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationCRN_pAF_BBLA.mat'], 's')
            
            case 'LAA'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationCRN_pAF_LAA';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationCRN_pAF_LAA.mat'], 's')

            case 'MVR'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationCRN_pAF_MVR';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationCRN_pAF_MVR.mat'], 's')
            
            case 'LA'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationCRN_pAF_LA';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationCRN_pAF_LA.mat'], 's')
        end

    case 'KOIVUMAKI'
        parfor k=1:trials
            [time{k},X{k},IsJs{k}] = simulationKOIV(p(k,:), 3, pAF_state);
        end
        
        % Save simulation and erase unnecessary elements
        switch atrial_region
            case 'PV'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationKOIV_pAF_PV';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationKOIV_pAF_PV.mat'], 's')
            case 'BBLA'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationKOIV_pAF_BBLA';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationKOIV_pAF_BBLA.mat'], 's')
            
            case 'LAA'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationKOIV_pAF_LAA';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationKOIV_pAF_LAA.mat'], 's')

            case 'MVR'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationKOIV_pAF_MVR';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationKOIV_pAF_MVR.mat'], 's')
            
            case 'LA'
                delete(poolobj);
                toc
                disp('Esperando a guardar');
                s.name='PopulationKOIV_pAF_LA';
                s.params=p;
                s.result=[(time'),(X'),(IsJs')];
                save(['PopulationKOIV_pAF_LA.mat'], 's')
        end
end
clearvars -except p trials
