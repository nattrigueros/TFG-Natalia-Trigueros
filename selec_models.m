%% SCRIPT FOR MODEL SELECTION INSIDE GENERATED POPULATION: 
% Outputs with the selected population matrix that fulfill all the 
% biomarkers physiological ranges.
% Natalia Trigueros
function [biomark,apd20_sum,apd50_sum,apd90_sum,apa_sum,dvdt_sum,rmp_sum,selec, sel1, sel2, sel3,s,ind_apd90,ind_rmp,ind_dvdtmax,ind_apa,ind_apd90_max] = selec_models(model_name)

% Selec is a matrix with a 1 for those model in which the biomarker
% fulfills the limits criteria and otherwise, a 0.
%
% sel1, sel2 and sel3 are the models that fulfill the callibration criteria
% If correctly executed, sel1 sel2 and sel3 must be EQUAL.
%
% Define the physiological range low and high limits for creating the POM
% Ranges are defined from experimental bibliography data
switch model_name
    case 'CRN CONTROL'
       f=load (['PopulationCRN_raw.mat']);
       s=f.s;
       lowlim=[0 95 184 88 152 -92];
       highlim=[15 183 360 119 288 -68];
    case 'CRN PV'
       f=load (['PopulationCRN_pAF_PV.mat']);
       s=f.s;
       lowlim=[0 112 160 88 152 -92];
       highlim=[15 179 285 119 288 -68];
    case 'CRN BBLA'
       f=load (['PopulationCRN_pAF_BBLA.mat']);
       s=f.s;
       lowlim=[0 112 160 88 152 -92];
       highlim=[15 179 285 119 288 -68];  
    case 'CRN LAA'
       f=load (['PopulationCRN_pAF_LAA.mat']);
       s=f.s;
       lowlim=[0 112 160 88 152 -92];
       highlim=[15 179 285 119 288 -68];
    case 'CRN MVR'
       f=load (['PopulationCRN_pAF_MVR.mat']);
       s=f.s;
       lowlim=[0 112 160 88 152 -92];
       highlim=[15 179 285 119 288 -68];
    case 'CRN LA'
       f=load (['PopulationCRN_pAF_LA.mat']);
       s=f.s;
       lowlim=[0 112 160 88 152 -92];
       highlim=[15 179 285 119 288 -68];  

    case 'KOIV CONTROL'
       f=load (['PopulationKOIV_raw.mat']);
       s=f.s;
       lowlim=[0 17.1 169.28 95.04 123.12 -86.48];
       highlim=[15.6 32.94 331.2 128.52 233.28 -63.92]; 
%        lowlim=[0 10 100 88 152 -92];
%        highlim=[15 183 360 140 288 -68];
    case 'KOIV PV'
       f=load (['PopulationKOIV_pAF_PV.mat']);
       s=f.s;
       lowlim=[0 20.16 89.6 101.2 139.84 -88.32];
       highlim=[14.7 32.22 159.6 136.85 264.96 -65.28]; 
    case 'KOIV BBLA'
       f=load (['PopulationKOIV_pAF_BBLA.mat']);
       s=f.s;
       lowlim=[0 20.16 89.6 101.2 139.84 -88.32];
       highlim=[14.7 32.22 159.6 136.85 264.96 -65.28];
    case 'KOIV LAA'
       f=load (['PopulationKOIV_pAF_LAA.mat']);
       s=f.s;
       lowlim=[0 20.16 89.6 101.2 139.84 -88.32];
       highlim=[14.7 32.22 159.6 136.85 264.96 -65.28];
    case 'KOIV MVR'
       f=load (['PopulationKOIV_pAF_MVR.mat']);
       s=f.s;
       lowlim=[0 20.16 89.6 101.2 139.84 -88.32];
       highlim=[14.7 32.22 159.6 136.85 264.96 -65.28];
    case 'KOIV LA'
       f=load (['PopulationKOIV_pAF_LA.mat']);
       s=f.s;
       lowlim=[0 20.16 89.6 101.2 139.84 -88.32];
       highlim=[14.7 32.22 159.6 136.85 264.96 -65.28];
end

% Get biomark vector
[biomark,apd20_sum,apd50_sum,apd90_sum,apa_sum,dvdt_sum,rmp_sum,s,ind_apd90,ind_rmp,ind_dvdtmax,ind_apa,ind_apd90_max] = calculateBiomarkers_calibration(s,model_name);

% Verify that all models fulfill the limits
selec= biomark(:,1:6,:) >lowlim & biomark(:,1:6,:)<highlim

% We can see which models are selected in each beat
% Beat 1
beat1 = selec(:,:,1);
suma=sum(beat1,2);
sel1=find(suma==6); % As mentioned, the 6 criteria must be fulfilled
% Beat 2
beat2 = selec(:,:,2);
suma=sum(beat2,2);
sel2=find(suma==6);
% Beat 3
beat3 = selec(:,:,3);
suma=sum(beat3,2);
sel3=find(suma==6);


% If for some mo del conditions are fulfilled in one beat, but not for other
% that beat is discarded (ALL 3 beats must fulfill conditions).
%
% Search models that have been selected in beat 1 but not in 2
ind=setdiff(sel1,sel2); 
for i=1:length(ind)
    % Eliminated from beat 1 the models not being selected in beat 2
    sel1(find(sel1==ind(i)))=[]; 
end
ind=setdiff(sel1,sel3);
for i=1:length(ind)
    sel1(find(sel1==ind(i)))=[];
end
ind=setdiff(sel2,sel1);
for i=1:length(ind)
    sel2(find(sel2==ind(i)))=[];
end
ind=setdiff(sel2,sel3);
for i=1:length(ind)
    sel2(find(sel2==ind(i)))=[];
end
ind=setdiff(sel3,sel1);
for i=1:length(ind)
    sel3(find(sel3==ind(i)))=[];
end
ind=setdiff(sel3,sel2);
for i=1:length(ind)
    sel3(find(sel3==ind(i)))=[];

end
end

