%% MAIN FOR GENERATING POPULATION OF PHYSIOLOGICAL MODELS
% CRN (Control/pAF), KOIVUMAKI(Control/pAF)
% Natalia Trigueros

%% STEP 1: Choose model from command window
% Choose the model selection thorugh the command window
m_n = input(['Choose model from options: a / b \n\na.CRN \nb.KOIVUMAKI \n\nYour selection is:'],'s');
reg = input(['Choose the kind of study you want to perform: a / b \n\na.CONTROL \nb.pAF in PV \nc.pAF in BBLA \nd.pAF in LAA \ne.pAF in MVR \nf.pAF in LA \n\nYour selection is:'],'s');

switch m_n
    % For CRN
    case {'a','a.','A','A.','CRN','a.CRN'}
        switch reg
            case {'a','A'}
                model_name='CRN CONTROL'
            case {'b','B'}
                model_name='CRN PV'
            case {'c','C'}
                model_name='CRN BBLA'
            case {'d','D'}
                model_name='CRN LAA'
            case {'e','E'}
                model_name='CRN MVR'
            case {'f','F'}
                model_name='CRN LA'
            otherwise 
                disp(['Run again file and choose a correct option.'])
        end

    % For Koivumaki            
    case {'b','b.','B','B.','KOIVUMAKI','b.KOIVUMAKI'}
        switch reg
            case {'a','A'}
                model_name='KOIV CONTROL'
            case {'b','B'}
                model_name='KOIV PV'
            case {'c','C'}
                model_name='KOIV BBLA'
            case {'d','D'}
                model_name='KOIV LAA'
            case {'e','E'}
                model_name='KOIV MVR'
            case {'f','F'}
                model_name='KOIV LA'
            otherwise 
                disp(['Run again file and choose a correct option a/b'])
        end
    
    otherwise 
        disp(['Run again file and choose a correct option a/b'])

end

%% STEP 2: Retrieve matrices of biomarker ranges fulfillment

% Get the matrix of 0/1s depicting models that fulfill conditions for 3
% last beats
[biomark,apd20_sum,apd50_sum,apd90_sum,apa_sum,dvdt_sum,rmp_sum,selec,sel1,sel2,sel3,s,ind_apd90,ind_rmp,ind_dvdtmax,ind_apa,ind_apd90_max] = selec_models(model_name);

% Show in command window the number of models that are in physiological
% range
disp(['We mantain ' num2str(size(sel3,1)) ' models for model ' num2str(model_name)])
disp([' '])

biomarkPhys1=zeros(length(sel1),6,1);
biomarkPhys2=zeros(length(sel2),6,1);
biomarkPhys3=zeros(length(sel3),6,1);

%% STEP 3: Save POM data in .mat
% Create matrices for storing the biomarkers values + currents/time/state
% vars vectors for each choosen as physiological model 
index = {};
s_r = {};
s_f = {};

for j=1:length(sel1)
    keep_j = sel1(j,1);
    index = [index;keep_j];
    biomark(keep_j,:,1);
    biomarkPhys1 = [biomarkPhys1;biomark(keep_j,:,1)];
end
b_p_1 = biomarkPhys1(length(sel1)+1:end,:);

for j=1:length(sel2)
    keep_j = sel1(j,1);
    biomark(keep_j,:,2);
    biomarkPhys2 = [biomarkPhys2;biomark(keep_j,:,2)];
end
b_p_2 = biomarkPhys2(length(sel2)+1:end,:);

for j=1:length(sel3)
    keep_j = sel3(j,1);
    biomark(keep_j,:,3);
    biomarkPhys3 = [biomarkPhys3;biomark(keep_j,:,3)];
end
b_p_3 = biomarkPhys3(length(sel3)+1:end,:);

for i=1:length(index)
    j = index{i};
    s_r = [s_r;s.result(j,:,:)];
end

for i=1:length(index)
    j = index{i};
    s_f = [s_f;s.params(j,:)];
end

% Store data in .mat
% Structure:
% - Params: Biomarkers for each beat. Columns 1-6 biomarkers third-too-last
% beat, columns 7-12 biomarker second-too-last beat,columns 13-18 
% biomarkers last beat.
% - Results: Rows = Models. First column time, second column state vars, 
% third column currents. Inside each three more columns corresponding to 
% each beat.

switch model_name
    case 'CRN CONTROL'
        sel1_CRN=sel1;
        disp('Esperando a guardar');
        t.name='PopPhysCRN_control';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysCRN_control.mat'], 't')
    case 'CRN PV'
        disp('Esperando a guardar');
        t.name='PopPhysCRN_pAF_PV';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysCRN_pAF_PV.mat'], 't')
    case 'CRN BBLA'
        disp('Esperando a guardar');
        t.name='PopPhysCRN_pAF_BBLA';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysCRN_pAF_BBLA.mat'], 't')
    case 'CRN LAA'
        disp('Esperando a guardar');
        t.name='PopPhysCRN_pAF_LAA';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysCRN_pAF_LAA.mat'], 't')
    case 'CRN MVR'
        disp('Esperando a guardar');
        t.name='PopPhysCRN_pAF_MVR';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysCRN_pAF_MVR.mat'], 't')
    case 'CRN LA'
        disp('Esperando a guardar');
        t.name='PopPhysCRN_pAF_LA';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysCRN_pAF_LA.mat'], 't')

    case 'KOIV CONTROL'
        sel1_KOIV=sel1;
        disp('Esperando a guardar');
        t.name='PopPhysKOIV_control';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysKOIV_control.mat'], 't')
    case 'KOIV PV'
        disp('Esperando a guardar');
        t.name='PopPhysKOIV_pAF_PV';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r; 
        save(['PopPhysKOIV_pAF_PV.mat'], 't')
    case 'KOIV BBLA'
        disp('Esperando a guardar');
        t.name='PopPhysKOIV_pAF_BBLA';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysKOIV_pAF_BBLA.mat'], 't')
    case 'KOIV LAA'
        disp('Esperando a guardar');
        t.name='PopPhysKOIV_pAF_LAA';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysKOIV_pAF_LAA.mat'], 't')
    case 'KOIV MVR'
        disp('Esperando a guardar');
        t.name='PopPhysKOIV_pAF_MVR';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysKOIV_pAF_MVR.mat'], 't')
    case 'KOIV LA'
        disp('Esperando a guardar');
        t.name='PopPhysKOIV_pAF_LA';
        t.params=s_f;
        t.biomarkers=[(b_p_1),(b_p_2),(b_p_3)];
        t.result=s_r;
        save(['PopPhysKOIV_pAF_LA.mat'], 't')
end

%% STEP 4: Calculate POM parameters
% Calculate mean value for all the studied biomarkers
%
%%%% Mean APD20
c_apd20_2 = 0;
for j=1:size(t.params,1) % Perform mean of 150 samples from each beat = 450 samples  
    c_apd20_1 = 0;
    for i=1:3
        apd20_1 = apd20_sum(j,1,i);
        c_apd20_1 = apd20_1 + c_apd20_1;
    end
    c_apd20_2 = c_apd20_2 + c_apd20_1;
end
av_apd20 = c_apd20_2/(size(t.params,1)*3);
disp(['The mean for APD20 is ' num2str(av_apd20) ' ms.'])

%%%% Mean APD50
c_apd50_2 = 0;
for j=1:size(t.params,1) % Perform mean of 150 samples from each beat = 450 samples    
    c_apd50_1 = 0;
    for i=1:3
        apd50_1 = apd50_sum(j,1,i);
        c_apd50_1 = apd50_1 + c_apd50_1;
    end
    c_apd50_2 = c_apd50_2 + c_apd50_1;
end
av_apd50 = c_apd50_2/(size(t.params,1)*3);
disp(['The mean for APD50 is ' num2str(av_apd50) ' ms.'])

%%%% Mean APD90
c_apd90_2 = 0;
for j=1:size(t.params,1) % Perform mean of 150 samples from each beat = 450 samples    
    c_apd90_1 = 0;
    for i=1:3
        apd90_1 = apd90_sum(j,1,i);
        c_apd90_1 = apd90_1 + c_apd90_1;
    end
    c_apd90_2 = c_apd90_2 + c_apd90_1;
end
av_apd90 = c_apd90_2/(size(t.params,1)*3);
disp(['The mean for APD90 is ' num2str(av_apd90) ' ms.'])

%%%% Mean APA
c_apa_2 = 0;
for j=1:size(t.params,1) % Perform mean of 150 samples from each beat = 450 samples  
    c_apa_1 = 0;
    for i=1:3
        apa_1 = apa_sum(j,1,i);
        c_apa_1 = apa_1 + c_apa_1;
    end
    c_apa_2 = c_apa_2 + c_apa_1;
end
av_apa = c_apa_2/(size(t.params,1)*3);
disp(['The mean for APA is ' num2str(av_apa) ' mV.'])

%%%% Mean dvdt_max
c_dvdt_2 = 0;
for j=1:size(t.params,1) % Perform mean of 150 samples from each beat = 450 samples   
    c_dvdt_1 = 0;
    for i=1:3
        dvdt_1 = dvdt_sum(j,1,i);
        c_dvdt_1 = dvdt_1 + c_dvdt_1;
    end
    c_dvdt_2 = c_dvdt_2 + c_dvdt_1;
end
av_dvdt = c_dvdt_2/(size(t.params,1)*3);
disp(['The mean for dvdt_max is ' num2str(av_dvdt) ' V/s.'])

%%%% Mean RMP
c_rmp_2 = 0;
for j=1:size(t.params,1) % Perform mean of 150 samples from each beat = 450 samples   
    c_rmp_1 = 0;
    for i=1:3
        rmp_1 = rmp_sum(j,1,i);
        c_rmp_1 = rmp_1 + c_rmp_1;
    end
    c_rmp_2 = c_rmp_2 + c_rmp_1;
end
av_rmp = c_rmp_2/(size(t.params,1)*3);
disp(['The mean for RMP is ' num2str(av_rmp) ' mV.'])

%% Save the indexes of the APs being discarded for resulting as Non-Physiological
switch model_name
    case 'CRN CONTROL'
        nonphys=zeros(1,1);
        pom=1:1000;
        log1= ismember(pom,sel1_CRN);
        
        for s=1:length(log1)
            memb1=log1(s);
            if memb1==0
                nonphys=[nonphys;s];
            end
        end
        nonphys=nonphys(2:end,:);
        
        nonphy.name='NonPhysCRN_AP';
        nonphy.index=nonphys;
        save(['NonPhysCRN_AP.mat'], 'nonphy')
        
        
        crit.name='NonPhysCRN_AP';
        crit.index=[ind_apd90,ind_rmp,ind_dvdtmax,ind_apa,ind_apd90_max];
        save(['critCRN_AP.mat'], 'crit')

    case 'KOIV CONTROL'
        nonphys=zeros(1,1);
        pom=1:1000;
        log2= ismember(pom,sel1_KOIV);
        
        for s=1:length(log2)
            memb2=log2(s);
            if memb2==0
                nonphys=[nonphys;s];
            end
        end
        nonphys=nonphys(2:end,:);
        
        nonphy.name='NonPhysKOIV_AP';
        nonphy.index=nonphys;
        save(['NonPhysKOIV_AP.mat'], 'nonphy')

        crit.name='NonPhysKOIV_AP';
        crit.index=[ind_apd90,ind_rmp,ind_dvdtmax,ind_apa,ind_apd90_max];
        save(['critKOIV_AP.mat'], 'crit')
end

