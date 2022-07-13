%% BIOMARKERS CALCULATION
% Outputs with biomarker vector for each trial.
% Natalia Trigueros

function [biomark,apd20_sum,apd50_sum,apd90_sum,apa_sum,dvdt_sum,rmp_sum,s,ind_apd90,ind_rmp,ind_dvdtmax,ind_apa,ind_apd90_max] = calculateBiomarkers_calibration(s,model_name)

% s is a simulation from the basal population.
%
% Every row from biomark is a member (trial) of the population (1000
% rows = 1000 trials).
%
% Each column in biomark is a biomarker. Biomarkers under study:
% - APD20  - APD50  - APD90  - APA  - dVdtmax  - RMP 
% In the following order: biomark=[apd20 apd90 apa dvdtmax rmp]
%             
% The third component makes reference to the beat (1=two before last; 2=
% one before last; 3=last)

% Create matrices for storing data
biomark=zeros(length(s.result),6,3);
apd20_sum=zeros(length(s.result),1,3);
apd50_sum=zeros(length(s.result),1,3);
apd90_sum=zeros(length(s.result),1,3);
apa_sum=zeros(length(s.result),1,3);
dvdt_sum=zeros(length(s.result),1,3);
rmp_sum=zeros(length(s.result),1,3);

val_apd90=100000;
val_apd90_max=0;
val_apa=100000;
val_dvdtmax=100000,
val_rmp=100000;

for j=1:length(s.result) % j goes from 1 to 1000 (nº of trials)
    
    for i=1:3 % corresponding to the number of beats under study
        
        % As models from CRN work with milliseconds and Koivumaki's with
        % seconds, two different functions have been created to treat 
        % population data. However, in the output (biomark matrix) BOTH 
        % models "time" data is in milliseconds.
        % Units:
        % biomark = [ms ms ms mV V/s mV]

        switch model_name
            case {'CRN CONTROL','CRN PV','CRN BBLA','CRN LAA','CRN MVR','CRN LA'}
                t=s.result{j,1}{1,i}; % Time 
                v=s.result{j,2}{1,i}(:,15); % Voltage
                I=s.result{j,3}{1,i}(:,[2:4,7:9]);

                % APD90
                [apd90,apa,dvdtmax,rmp,vpeak]=APD_CRN(t,v,90);
                % APD20
                apd20=APD_CRN(t,v,20);
                % APD50
                apd50=APD_CRN(t,v,50);
                % Construct biomark vector for every member
                %                1      2      3   4     5     6   
                biomark(j,:,i)=[apd20 apd50 apd90 apa dvdtmax rmp];
                if i==1
                    if apd90<val_apd90
                        val_apd90=apd90;
                        ind_apd90=j;
                    end
                    if apd90>val_apd90_max
                        val_apd90_max=apd90;
                        ind_apd90_max=j;
                    end
                    if apa<val_apa
                        val_apa=apa;
                        ind_apa=j;
                    end
                    if dvdtmax<val_dvdtmax
                        val_dvdtmax=dvdtmax;
                        ind_dvdtmax=j;
                    end
                    if rmp<val_rmp
                        val_rmp=rmp;
                        ind_rmp=j;   
                    end
                end 
                
            case {'KOIV CONTROL','KOIV PV','KOIV BBLA','KOIV LAA','KOIV MVR','KOIV LA'}
                t=s.result{j,1}{1,i};  % Time 
                v=s.result{j,2}{1,i}(:,1); % Voltage
                I=s.result{j,3}{1,i}(:,[2:4,7:9]);

                % APD90
                [apd90,apa,dvdtmax,rmp,vpeak]=APD_KOIV(t,v,90);
                % APD20
                apd20=APD_KOIV(t,v,20);
                % APD50
                apd50=APD_KOIV(t,v,50);
                % Construct biomark vector for every member
                %                1      2      3   4     5     6   
                biomark(j,:,i)=[apd20 apd50 apd90 apa dvdtmax rmp];
                if i==1
                    if apd90<val_apd90
                        val_apd90=apd90;
                        ind_apd90=j;
                    end
                    if apd90>val_apd90_max
                        val_apd90_max=apd90;
                        ind_apd90_max=j;
                    end
                    if apa<val_apa
                        val_apa=apa;
                        ind_apa=j;
                    end
                    if dvdtmax<val_dvdtmax
                        val_dvdtmax=dvdtmax;
                        ind_dvdtmax=j;
                    end
                    if rmp<val_rmp
                        val_rmp=rmp;
                        ind_rmp=j;   
                    end
                end 
        end
        
    end
 
end

% This is for calculations purposes of the average values for each
% biomarker.
apd20_sum=apd20_sum+biomark(:,1,:);
apd50_sum=apd50_sum+biomark(:,2,:);
apd90_sum=apd90_sum+biomark(:,3,:);
apa_sum=apa_sum+biomark(:,4,:);
dvdt_sum=dvdt_sum+biomark(:,5,:);
rmp_sum=rmp_sum+biomark(:,6,:);

end
%% Calculation of biomarkers function: COURTEMANCHE
function [APDxx,Vpp,dVdt_max,Vr,Vp]=APD_CRN(t,v,xx)

%APD at different repolarization
%Other parameters: Amplitude (Vpp), max upstroke (dVdt_max), Vrest, Vpeak
    
[Vp, Vp_index]=max(v); %Caution with elevated domes
% Get RMP value:
[Vr, Vr_index]=min(v(1:end));
% Get APA value:
Vpp=Vp-Vr;

% Derivative of AP
dVdt=zeros(length(v)-1,1);
for i=1:(length(v)-1)
    dVdt(i)=(v(i+1)-v(i))/(t(i+1)-t(i));
end
% Detect EADs
t_ind_repol= find(t>150,1); % Time index 150 after initial time
if sum(dVdt(t_ind_repol:end)>0.05)>1  
    if Vr > -70
        APDxx=NaN; Vpp=NaN; dVdt_max=NaN; Vr=NaN; Vp=NaN;
    end
else
    % Get dVdtmax value:
    [dVdt_max, dVdt_max_index]=max(dVdt);
    t_maxderv=t(dVdt_max_index);

    V_xxrep=Vp-(xx/100)*Vpp;
    v_rep=v(Vp_index:end); % Attention: Now Vp is index 1. I fix it here (*)
    [dif,vxx_index_cut]=min(abs(v_rep-V_xxrep));
    vxx_index=vxx_index_cut+(Vp_index-1); %(*)

    %Interpolation in order to get close to V_xxrep (it isn't a real point):
    t1=t(vxx_index); v1=v(vxx_index);
    if v(vxx_index)>V_xxrep
        t2=t(vxx_index+1); v2=v(vxx_index+1);
    end
    if v(vxx_index)<V_xxrep
        t2=t(vxx_index-1); v2=v(vxx_index-1);
    end
    t_vxx=t1+(V_xxrep-v1)*((t2-t1)/(v2-v1));

    % Get final APD value:
    APDxx=t_vxx-t_maxderv;
end
end

%% Calculation of biomarkers function: KOIVUMAKI
function [APDxx,Vpp,dVdt_max,Vr,Vp]=APD_KOIV(t,v,xx)

%APD at different repolarization
%Other parameters: Amplitude (Vpp), max upstroke (dVdt_max), Vrest, Vpeak
    
[Vp, Vp_index]=max(v); % Caution with elevated domes
% Get RMP value:
[Vr, Vr_index]=min(v(1:end));
% Get APA value:
Vpp=Vp-Vr;

%Derivative of AP
dVdt=zeros(length(v)-1,1);
for i=1:(length(v)-1)
    dVdt(i)=(v(i+1)-v(i))/(1000*t(i+1)-1000*t(i));
end

% Get dVdtmax value:
[dVdt_max, dVdt_max_index]=max(dVdt);
t_maxderv=t(dVdt_max_index);

V_xxrep=Vp-(xx/100)*Vpp;
v_rep=v(Vp_index:end); % Attention: Now Vp is index 1. I fix it here (*)
[dif,vxx_index_cut]=min(abs(v_rep-V_xxrep));
vxx_index=vxx_index_cut+(Vp_index-1); % (*)

% Interpolation in order to get close to V_xxrep (it isn't a real point):
t1=t(vxx_index); v1=v(vxx_index);
if v(vxx_index)>V_xxrep
    t2=t(vxx_index+1); v2=v(vxx_index+1);
end
if v(vxx_index)<V_xxrep
    t2=t(vxx_index-1); v2=v(vxx_index-1);
end
t_vxx=t1*1000+(V_xxrep-v1)*((t2-t1)*1000/(v2-v1));

% Get final APD value:
APDxx=t_vxx-t_maxderv;
end