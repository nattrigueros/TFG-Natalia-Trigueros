%% KOIVUMAKI MODEL 2011
% With added pAF remodelling + IACh current expression
% Natalia Trigueros

function [dy, currents, t] = Koivumaki_model(t, y,ACh, mutation, scaling_factors, pAF_state)
%% ************************************************************************
%  Human atrial myocyte model
%
%  JT Koivumaeki, T Korhonen, P Tavi (2011)
%  Impact of Sarcoplasmic Reticulum Calcium Release on Calcium Dynamics and
%  Action Potential Morphology in Human Atrial Myocytes: A Computational
%  Study.
%  PLoS Comput Biol 7(1): e1001067.

%  http://dx.doi.org/10.1371/journal.pcbi.1001067
%
%  - Model developed based on:
%    * Nygren et al. (1998) model, including
%**************************************************************************

%% Additional Parameters

% pAF remodelling
r_kr = 1;
r_CaL = 1; 
r_to = 1;
r_ks = 1;
r_K1 = 1;
switch pAF_state

    case 'pAF_no'   % Control
        rik1=1;
        rikACh=1;

    case 'pAF_yes'  % pAF remodelling
        rik1=1.5;
        rikACh=1.5;
 
end

switch ACh   % ACh concentration: Taking 0.005 as baseline concentration.
    case "ACh005"
        ACh = 0.005;
    case "ACh0"
        ACh = 0.0;
end

switch mutation
    case "NO"
        p0r  =  1;
        p1r  =  1;
        p2r  =  1;
        p5r  =  1;
        p6r  =  1;
        p7r  =  1;
        p0t  =  1;
        p5t  =  1;
    case "WT"
        p0r  =  1;
        p1r  =  0;
        p2r  =  1;
        p5r  =  0;
        p6r  =  1;
        p7r  =  1;
        p0t  =  0;
        p5t  =  1;
        
    case "T895M"
        
        p0r  =  1.405;
        p1r  =  0;
        p2r  =  1.02;
        p5r  =  0;
        p6r  =  1;
        p7r  =  2.043;
        p0t  =  0;
        p5t  =  1;
        
    case "T436M"
        p0r  =  1.513;
        p1r  =  0;
        p2r  =  1.13;
        p5r  =  0;
        p6r  =  1;
        p7r  =  1.394;
        p0t  =  0;
        p5t  =  1;
        
    case "V17M"
        p0r  =  1.956;
        p1r  =  0;
        p2r  =  1;
        p5r  =  1.786;
        p6r  =  .973;
        p7r  =  5;
        p0t  =  13.04;
        p5t  =  7.587663;
        
end

%% Definition of differential variables
i_V = 1;
i_lastend = i_V;

% INa
i_start = i_lastend + 1;
i_INam = i_start;       i_INah1 = i_start + 1;      i_INah2 = i_start + 2;
i_lastend = i_INah2;

% ICaL 
i_start = i_lastend + 1;
i_ICaLd = i_start;      i_ICaLf1 = i_start + 1;         i_ICaLf2 = i_start + 2;      i_ICaLfca = i_start + 3;
i_lastend = i_ICaLfca;

% It
i_start = i_lastend + 1;
i_Itr = i_start;        i_Its = i_start + 1;
i_lastend = i_Its;

% Isus (Ikur)
i_start = i_lastend + 1;
i_Isusr = i_start;       i_Isuss = i_start + 1;
i_lastend = i_Isuss;

% IKs
i_start = i_lastend + 1;
i_IKsn = i_start;
i_lastend = i_IKsn;

% IKr
i_start = i_lastend + 1;
i_IKrpa = i_start;
i_lastend = i_IKrpa;

% If
i_start = i_lastend + 1;
i_Ify = i_start;
i_lastend = i_Ify;

% RyR
i_start = i_lastend + 1;
i_RyRoss = i_start; i_RyRcss = i_start + 1; i_RyRass = i_start + 2;
i_RyRo1 = i_start + 3; i_RyRc1 = i_start + 4; i_RyRa1 =  i_start + 5;
i_RyRo2 = i_start + 6; i_RyRc2 = i_start + 7; i_RyRa2 =  i_start + 8;
i_RyRo3 = i_start + 9; i_RyRc3 = i_start + 10; i_RyRa3 =  i_start + 11;
i_lastend = i_RyRa3;

% SERCA
i_start = i_lastend + 1;
i_SERCACa1 = i_start; % first compartment from center
i_SERCACa2 = i_start + 1; % 2nd compartment from center
i_SERCACa3 = i_start + 2; % 3rd compartment from center
i_SERCACass = i_start + 3; % ss-SERCA
i_lastend = i_SERCACass;

% Nai ja Ki
i_start = i_lastend + 1;
i_Nai = i_start;
i_Ki = i_start + 1;
i_lastend = i_Ki;

% Cai and CaSR 
i_start = i_lastend + 1;
i_Cass = i_start;
i_Cacenter = i_start + 1;
%i_lastend = i_Cacenter;


%% Select stimulus ('stim_I_adjustable_duration_and_amplitude')
stim_offset = 0.0; % has to be smaller than the max step in solver
stim_duration = 0.001; % in seconds
stim_amplitude = -2200;% (44pA/pF)
stim_BCL = 1000;
stim_period = stim_BCL/1000;
if ((t-floor(t/stim_period)*stim_period >= stim_offset) && (t-floor(t/stim_period)*stim_period <= stim_offset + stim_duration))
    Istim = scaling_factors(2)*stim_amplitude;
else
    Istim = 0.0;
end
        
%% Physical & environmental constants 
F = 96487; %(Table 1 Koivumaki et al. 2011)
R = 8314; %(Table 1 Koivumaki et al. 2011)
T = 33 + 273.15; % 306.15 K (Table 1 Koivumaki et al. 2011)
Cm = 0.05; % nF (Table 1 Koivumaki et al. 2011)

% Geometry (Table 1 Koivumaki et al. 2011)
Vss = 4.99232e-5; % nL 
rjunct = 6.5; % mum 
lcell = 122.051; % mum
deltar = 1.625; % mum

% Ca diffusion grid (Text S1 Koivumaki et al. 2011)
Aj_nj = pi*rjunct*2*lcell*0.5; % Area between junct and nonjunct (eq.16)
xj_nj = 0.02/2 + deltar/2; % diffusion distance from center to center of junct to first njunct (eq. 17)

% Diffusion compartment volumes (Text S1 Koivumaki et al. 2011)
j=1:4;
Vnonjunct = zeros(length(j),1);
Vnonjunct = (pi.*(j.*deltar).^2.*lcell-pi.*((j-1).*deltar).^2.*lcell).*1e-6.*0.5; %nL (eq. 42)
VSR = 0.0225.*Vnonjunct; %(eq. 43)
Vcytosol = sum(Vnonjunct) + Vss;

% Non-junct Cai data & CaSR data
Cai = zeros(length(j),1);
Cai = y(i_Cacenter:length(j)+i_Cacenter-1);
CaSR = y(length(j)+i_Cacenter:length(j)*2+i_Cacenter-1);

% Cytosol Ca Buffers (Table 2 Koivumaki et al. 2011)
BCa = 0.024; %mM 
SLlow = 165; %mM
SLhigh = 13; %mM
KdBCa = 0.00238; %mM
KdSLlow = 1.1; %mM
KdSLhigh = 0.013; %mM

% SR Ca buffers
CSQN =  6.7;
KdCSQN = 0.8;


%% Ion channel conductances & permeabilities & other parameters
gNa = scaling_factors(3) * 0.00182; %(Table S1 Koivumaki et al. 2011)
gCaL = scaling_factors(4) * 25.3125 * r_CaL; % (p-7 Text S1 Koivumaki et al. 2011)
gt = scaling_factors(5) * 7.5 * r_to; %(Table 18 Nyegren et al. 1998) nS
gsus = scaling_factors(6) *2.75; %(Table 18 Nyegren et al. 1998)
gKs = scaling_factors(8) *1.0 * r_ks; %(Table 18 Nyegren et al. 1998)
gKr = scaling_factors(7) *0.5 * r_kr; %(Table 18 Nyegren et al. 1998) nS
gK1 = scaling_factors(9) *3.825 * r_K1; %(Table S1 Koivumaki et al. 2011)
gNab = scaling_factors(11) *0.060599; %(Table 18 Nyegren et al. 1998)
gCab = scaling_factors(12) *0.0952; %(Table S1 Koivumaki et al. 2011)
gIf = scaling_factors(10)*1; % (Zorn-Pauly 2004 LAW fit)

ECa_app=60; %mV %(Table 17 Nyegren et al. 1998)
INaKmax = scaling_factors(14)*70.8253; %(Table 17 Nyegren et al. 1998)
kNaKK = 1; %(Table 17 Nyegren et al. 1998)
kNaKNa = 11; %(Table 17 Nyegren et al. 1998)
ICaPmax = 2.0; %(Table S1 Koivumaki et al. 2011)
kCaP = scaling_factors(13) *0.0002; %(Table S1 Koivumaki et al. 2011)
kNaCa = scaling_factors(15)*0.008433945; %(Table S1 Koivumaki et al. 2011 (gINCX))
gam = 0.45; %(Table 17 Nyegren et al. 1998)
dNaCa = 0.0003; %(Table 17 Nyegren et al. 1998)

% Ca diffusion (Table 2 Koivumaki et al. 2011)
DCa = 780; 
DCaSR = 44; 
DCaBm = 25; 

% SERCA
SERCAKmf = 0.00025; %mM (Table S1 Koivumaki et al. 2011)
SERCAKmr = 1.8; %mM (Table S1 Koivumaki et al. 2011)

k4 = 7.5; % (Table S1 Koivumaki et al. 2011)
k3 = k4/SERCAKmr^2; % (eq 3 Text S1 Koivumaki et al. 2011)
k1 = 10^6 * k4; % (eq 1 Text S1 Koivumaki et al. 2011)
k2 = k1 * SERCAKmf^2; %(eq 2 Text S1 Koivumaki et al. 2011)
cpumps = 0.04; %mM (Table S1 Koivumaki et al. 2011)

% RyR (Table S2 Koivumaki et al. 2011)
RyRtauadapt = 1;
RyRtauactss = 5e-3; 
RyRtauinactss = 15e-3; 
RyRtauact = 18.75e-3; 
RyRtauinact = 87.5e-3;

% SR Ca leak
kSRleak = 0.006; %(Table S1 Koivumaki et al. 2011)

% Concentrations
Nao = 130; Ko = 5.4; Cao = 1.8; % Table 1 Koivumaki et al. 2011
Nai = y(i_Nai); Ki = y(i_Ki); Cass = y(i_Cass);

%% Analytical equations

% Reversal potentials
ENa = R*T/F * log ( Nao / Nai );
EK = R*T/F * log ( Ko / Ki );
ECa = R*T/F/2 * log ( Cao / Cass );

% INa (Nygren et al. 1998)
PNa = gNa;
INa = PNa .* y(i_INam).^3 .* ( 0.9.*y(i_INah1) + 0.1.*y(i_INah2) ) * Nao * y(i_V) * F^2/(R*T) * ( exp( (y(i_V)-ENa)*F/R/T ) - 1) / ( exp( y(i_V)*F/R/T ) - 1);
INaminf = 1/(1+exp((y(i_V)+27.12)/-8.21));
INahinf = 1/(1+exp((y(i_V)+63.6)/5.3));
INah1inf = INahinf;
INah2inf = INahinf;
INamtau = 0.000042*exp( -((y(i_V)+25.57)/28.8).^2 ) + 0.000024;
INah1tau = 0.03/(1+exp((y(i_V)+35.1)/3.2)) + 0.0003;
INah2tau = 0.12/(1+exp((y(i_V)+35.1)/3.2)) + 0.003;

% ICaL
ICaLfcainf = 1-(1/(1+(0.001/Cass)^2)); %(eq 45 Text S1 Koivumaki et al. 2011)
ICaL = gCaL * y(i_ICaLd) * y(i_ICaLf1) * y(i_ICaLf2) * y(i_ICaLfca) * (y(i_V) - ECa_app); %(eq 44 Text S1 Koivumaki et al. 2011)
ICaLdinf = 1/(1+exp((y(i_V)+9)/-5.8)); % Nygren et al.
ICaLfinf = 1/(1+exp((y(i_V)+27.4)/7.1)); % Nygren et al.
ICaLdtau = (0.0027 * exp(-((y(i_V) + 35)/30).^2) + 0.002); % Nygren et al.
ICaLf1tau = 0.98698.*exp(-((y(i_V)+30.16047)./7.09396).^2 ) + 0.04275./(1+exp((y(i_V)-51.61555)/-80.61331)) + 0.03576./(1+exp((y(i_V)+29.57272)/13.21758))-0.00821; % PLoS 2011
ICaLf2tau = 1.3323*exp(-((y(i_V)+40.0)/14.2)^2.0)+0.0626; % Nygren et al.
ICaLfcatau = 2e-3; % Ten Tusscher 

% Ito (Nygren et al. 1998)
It = gt * p5t * y( i_Itr ) * y( i_Its ) * ( y(i_V) - EK );
Itrinf = 1 / ( 1 + exp( ( y( i_V ) - 1) / -11) ); 
Itsinf = 1 / ( 1 + exp( ( y( i_V ) + 40.5 + p0t ) / 11.5) ); 
Itrtau = 0.0035 * exp( -( ( y( i_V ) / 30).^ 2 ) ) + 0.0015; 
Itstau = 0.4812 * exp( -( ( y( i_V ) + 52.45 ) / 14.97).^ 2 ) + 0.01414; 

% Isus / IKur (Nygren et al. (1998))
Isus = gsus * y(i_Isusr) * y(i_Isuss) * (y(i_V) - EK); 
Isusrinf = 1/(1 + exp((y(i_V) + 4.3)/-8.0)); 
Isussinf = 0.4/(1 + exp((y(i_V) + 20.0)/10)) + 0.6; 
Isusrtau = 0.009/(1 + exp((y(i_V) + 5)/12)) + 0.0005; 
Isusstau = 0.047/(1 + exp((y(i_V) + 60)/10)) + 0.3; 

% IKr  (Nygren et al. (1998))
IKrpi = 1 / ( 1 + exp( ( y( i_V ) + 55 + p5r ) / (24 * p6r) )); 
IKr = p7r * gKr * y(i_IKrpa) * IKrpi * ( y( i_V ) - EK ); 
IKrpainf = 1 / ( 1 + exp( ( y( i_V ) + 15 ) / -6) );                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
IKrpatau = p0r * 0.21718 * exp( -( ( y( i_V ) + 20.1376 + p1r ) / (22.1996 * p2r) ).^ 2 ) + 0.03118;

% IKs (Nygren et al. (1998))
IKs = gKs * y(i_IKsn) * (y(i_V) - EK); 
IKsninf = 1/(1+exp((y(i_V)-19.9)/-12.7));
IKsntau = 0.4*exp( -((y(i_V)-20)/20).^2 ) + 0.7;

% IK1 (Nygren et al. (1998))
IK1 = rik1 * gK1 * Ko.^0.4457 * (y(i_V) - EK) / (1 + exp(1.5*(y(i_V)-EK+3.6)*F/R/T));

% IKACh from Grandi - Rebecca 05.07.2021
i_KACh = rikACh * 1/(1+(0.03/ACh)^(2.1)) * (y(i_V) - EK) * (0.08 + 0.0400 / (1 + exp((y(i_V) + 91)/12)));

% Background leaks (Nygren et al. (1998))
INab = gNab * (y(i_V) - ENa);
ICab = gCab * (y(i_V) - ECa);

% INaK (Nygren et al. (1998))
INaK = INaKmax * Ko./(Ko + kNaKK) * Nai.^1.5/(Nai.^1.5 + kNaKNa.^1.5) * (y(i_V) + 150) / (y(i_V) + 200);

% INaCa (Nygren et al. (1998))
INaCa = kNaCa * ( (exp( gam.*y(i_V)*F/R/T ) .* Nai.^3 .* Cao - exp( (gam-1).*y(i_V)*F/R/T ) .* Nao^3 .* Cass) / ( 1 + dNaCa*(Nao^3 .* Cass + Nai.^3 .* Cao) ) );

% ICaP (Nygren et al. (1998))
ICaP = ICaPmax * Cass / (kCaP + Cass);

% If, Zorn-Pauly LAW fit
Ifyinf = 1 / (1 + exp((y(i_V)+97.82874)/12.48025));
Ifytau = 1 ./ (0.00332.*exp(-y(i_V)./16.54103)+23.71839.*exp(y(i_V)./16.54103));
IfNa = gIf * y(i_Ify)*((0.2677)*(y(i_V)-ENa));
IfK = gIf * y(i_Ify)*((1-0.2677)*(y(i_V)-EK)); 
If = IfK + IfNa;

% Ca buffers (Text S1 eq 28-31 Koivumaki et al 2011)
betass = ( 1 + SLlow*KdSLlow./(Cass + KdSLlow).^2 + SLhigh*KdSLhigh./(Cass + KdSLhigh).^2 + BCa*KdBCa./(Cass + KdBCa).^2  ).^(-1);
betai = ( 1 + BCa.*KdBCa./(Cai + KdBCa).^2  ).^(-1);
gammai = BCa.*KdBCa./(Cai + KdBCa).^2;
betaSR = ( 1 + CSQN.*KdCSQN./(CaSR + KdCSQN).^2 ).^(-1); 

%  Diffusion from junct to non-junct  (Text S1 eq 18 Koivumaki et al 2011)
Jj_nj = DCa * Aj_nj / xj_nj * (Cass-Cai(end)).*1e-6;

% SERCA fluxes  (Text S1 eq 4-5 Koivumaki et al 2011)
J_SERCASR1 = (-k3*CaSR(1).^2*(cpumps-y(i_SERCACa1))+k4*y(i_SERCACa1))*Vnonjunct(1)*2; % in 1 nl volume
J_bulkSERCA1 = (k1*Cai(1).^2*(cpumps-y(i_SERCACa1))-k2*y(i_SERCACa1))*Vnonjunct(1)*2; % in 1 nl volume

J_SERCASR2 = (-k3*CaSR(2).^2*(cpumps-y(i_SERCACa2))+k4*y(i_SERCACa2))*Vnonjunct(2)*2; % in 1 nl volume
J_bulkSERCA2 = (k1*Cai(2).^2*(cpumps-y(i_SERCACa2))-k2*y(i_SERCACa2))*Vnonjunct(2)*2; % in 1 nl volume

J_SERCASR3 = (-k3*CaSR(3).^2*(cpumps-y(i_SERCACa3))+k4*y(i_SERCACa3))*Vnonjunct(3)*2; % in 1 nl volume
J_bulkSERCA3 = (k1*Cai(3).^2*(cpumps-y(i_SERCACa3))-k2*y(i_SERCACa3))*Vnonjunct(3)*2; % in 1 nl volume

J_SERCASRss = (-k3*CaSR(4).^2*(cpumps-y(i_SERCACass))+k4*y(i_SERCACass))*Vss*2; % in 1 nl volume
J_bulkSERCAss = (k1*Cass.^2*(cpumps-y(i_SERCACass))-k2*y(i_SERCACass))*Vss*2; % in 1 nl volume

% RyR (Text S1 eq 7-11 Koivumaki et al 2011)
nuss = 625*Vss; %Table S2 Koivumaki et al. 2011
RyRSRCass = (1 - 1./(1 +  exp((CaSR(4)-0.3)./0.1)));
RyRainfss = 0.505-0.427./(1 + exp((Cass.*1000-0.29)./0.082));
RyRoinfss = (1 - 1./(1 +  exp((Cass.*1000-(y(i_RyRass) + 0.22))./0.03)));
RyRcinfss = (1./(1 + exp((Cass.*1000-(y(i_RyRass)+0.02))./0.01)));
Jrelss = nuss * ( y(i_RyRoss) ) * y(i_RyRcss) * RyRSRCass * ( CaSR(4) -  Cass ); 

nu1 = 1*Vnonjunct(1); %Table S2 Koivumaki et al. 2011
RyRSRCa1 = (1 - 1./(1 +  exp((CaSR(1)-0.3)./0.1)));
RyRainf1 = 0.505-0.427./(1 + exp((Cai(1).*1000-0.29)./0.082));
RyRoinf1 = (1 - 1./(1 +  exp(( Cai(1).*1000-(y(i_RyRa1) + 0.22))./0.03)));
RyRcinf1 = (1./(1 +  exp(( Cai(1).*1000-(y(i_RyRa1)+0.02))./0.01)));
Jrel1 = nu1 * ( y(i_RyRo1) ) * y(i_RyRc1) * RyRSRCa1 * ( CaSR(1) -  Cai(1) ); 

nu2 = 1*Vnonjunct(2); %Table S2 Koivumaki et al. 2011
RyRSRCa2 = (1 - 1./(1 +  exp((CaSR(2)-0.3)./0.1)));
RyRainf2 =  0.505-0.427./(1 + exp((Cai(2).*1000-0.29)./0.082));
RyRoinf2 = (1 - 1./(1 +  exp(( Cai(2).*1000-(y(i_RyRa2) + 0.22))./0.03)));
RyRcinf2 = (1./(1 +  exp(( Cai(2).*1000-(y(i_RyRa2)+0.02))./0.01)));
Jrel2 = nu2 * ( y(i_RyRo2) ) * y(i_RyRc2) * RyRSRCa2 * ( CaSR(2) -  Cai(2) ); 

nu3 = 1*Vnonjunct(3); %Table S2 Koivumaki et al. 2011
RyRSRCa3 = (1 - 1./(1 +  exp((CaSR(3)-0.3)./0.1)));
RyRainf3 =  0.505-0.427./(1 + exp((Cai(3).*1000-0.29)./0.082));
RyRoinf3 = (1 - 1./(1 +  exp(( Cai(3).*1000-(y(i_RyRa3) + 0.22))./0.03)));
RyRcinf3 = (1./(1 +  exp(( Cai(3).*1000-(y(i_RyRa3)+0.02))./0.01)));
Jrel3 = nu3 * ( y(i_RyRo3) ) * y(i_RyRc3) * RyRSRCa3 * ( CaSR(3) -  Cai(3) ); 

% SR leak fluxes (Text S1 eq 15 Koivumaki et al 2011)
JSRCaleak1 = kSRleak * ( CaSR(1) - Cai(1) ) * Vnonjunct(1); 
JSRCaleak2 = kSRleak * ( CaSR(2) - Cai(2) ) * Vnonjunct(2); 
JSRCaleak3 = kSRleak * ( CaSR(3) - Cai(3) ) * Vnonjunct(3); 
JSRCaleakss = kSRleak * ( CaSR(4) - Cass ) * Vss;

% Cafluxes in 1 nl volume (Text S1 eq 19-27 Koivumaki et al 2011)
JCa = zeros(length(j),1);
JCa(1) = -J_bulkSERCA1 + JSRCaleak1 + Jrel1;
JCa(2) = -J_bulkSERCA2 + JSRCaleak2 + Jrel2;
JCa(3) = -J_bulkSERCA3 + JSRCaleak3 + Jrel3;
JCa(4) = Jj_nj;
JCass = -Jj_nj + JSRCaleakss - J_bulkSERCAss + Jrelss;

JSRCa = zeros(length(j),1);
JSRCa(1) = J_SERCASR1 - JSRCaleak1 - Jrel1;
JSRCa(2) = J_SERCASR2 - JSRCaleak2 - Jrel2;
JSRCa(3) = J_SERCASR3 - JSRCaleak3 - Jrel3;
JSRCa(4) = J_SERCASRss - JSRCaleakss - Jrelss;


%% Differential equations
dy = zeros(length(j)*2+i_Cacenter-1,1);

% V
Iion= scaling_factors(1)*(INa + ICaL + It + Isus + IK1  + IKr + IKs + INab + ICab + INaK + ICaP + INaCa + If + i_KACh );
dy(i_V) = (Iion + Istim)/(-Cm); % currents are in (pA):

% INa
dy(i_INam) = (INaminf - y(i_INam))/INamtau;
dy(i_INah1) = (INah1inf - y(i_INah1))/INah1tau;
dy(i_INah2) = (INah2inf - y(i_INah2))/INah2tau;

% ICaL
dy(i_ICaLd) = (ICaLdinf - y(i_ICaLd))/ICaLdtau;
dy(i_ICaLf1) = (ICaLfinf - y(i_ICaLf1))/ICaLf1tau;
dy(i_ICaLf2) = (ICaLfinf - y(i_ICaLf2))/ICaLf2tau;
dy(i_ICaLfca) = (ICaLfcainf - y(i_ICaLfca))/ICaLfcatau;

% It
dy(i_Itr) = (Itrinf - y(i_Itr))/Itrtau;
dy(i_Its) = (Itsinf - y(i_Its))/Itstau;

% Isus
dy(i_Isusr) = (Isusrinf - y(i_Isusr))/Isusrtau;
dy(i_Isuss) = (Isussinf - y(i_Isuss))/Isusstau;

% IKs
dy(i_IKsn) = (IKsninf - y(i_IKsn))/IKsntau;

% IKr
dy(i_IKrpa) = (IKrpainf - y(i_IKrpa))/IKrpatau;

% If
dy(i_Ify) = (Ifyinf - y(i_Ify))/Ifytau;


% SERCACa
dy(i_SERCACa1) = 0.5*(-J_SERCASR1 + J_bulkSERCA1)/Vnonjunct(1); 
dy(i_SERCACa2) = 0.5*(-J_SERCASR2 + J_bulkSERCA2)/Vnonjunct(2); 
dy(i_SERCACa3) = 0.5*(-J_SERCASR3 + J_bulkSERCA3)/Vnonjunct(3); 
dy(i_SERCACass) = 0.5*(-J_SERCASRss + J_bulkSERCAss)/Vss; 

% RyR
dy(i_RyRoss) = (RyRoinfss-y(i_RyRoss))./RyRtauactss;
dy(i_RyRcss) = (RyRcinfss-y(i_RyRcss))./RyRtauinactss;
dy(i_RyRass) = (RyRainfss-y(i_RyRass))./RyRtauadapt;
dy(i_RyRo1) = (RyRoinf1-y(i_RyRo1))./RyRtauact;
dy(i_RyRc1) = (RyRcinf1-y(i_RyRc1))./RyRtauinact;
dy(i_RyRa1) = (RyRainf1-y(i_RyRa1))./RyRtauadapt;
dy(i_RyRo2) = (RyRoinf2-y(i_RyRo2))./RyRtauact;
dy(i_RyRc2) = (RyRcinf2-y(i_RyRc2))./RyRtauinact;
dy(i_RyRa2) = (RyRainf2-y(i_RyRa2))./RyRtauadapt;
dy(i_RyRo3) = (RyRoinf3-y(i_RyRo3))./RyRtauact;
dy(i_RyRc3) = (RyRcinf3-y(i_RyRc3))./RyRtauinact;
dy(i_RyRa3) = (RyRainf3-y(i_RyRa3))./RyRtauadapt;

% Nai & Ki
dy(i_Nai) = -(INa + INab + 3*INaK + 3*INaCa + IfNa)/ (Vcytosol*F); % currents are in pA
dy(i_Ki) = -(It + Isus + IK1 + IKr + IKs - 2*INaK + IfK + Istim + i_KACh) / (Vcytosol*F); % currents are in pA

% Ca
dy(i_Cass) = betass * ( JCass/Vss + (-ICaL - ICab - ICaP + 2*INaCa) / (2*Vss*F) ); % currents are in pA

dx=deltar;
dCaidt = zeros(length(Cai),1); 
dCaidt(1) = betai(1) .* (DCa + gammai(1).*DCaBm) .* ( (Cai(2)-2.*Cai(1)+Cai(1))./dx.^2 + (Cai(2)-Cai(1))./(2.*j(1).*dx.^2) ) - 2.*betai(1).*gammai(1).*DCaBm./(KdBCa + Cai(1)) .* ((Cai(2)-Cai(1))./(2.*dx)).^2 + JCa(1)./Vnonjunct(1).*betai(1); 
dCaidt(2) = betai(2) .* (DCa + gammai(2).*DCaBm) .* ( (Cai(3)-2.*Cai(2)+Cai(1))./dx.^2 + (Cai(3)-Cai(1))./(2.*j(2).*dx.^2) ) - 2.*betai(2).*gammai(2).*DCaBm./(KdBCa + Cai(2)) .* ((Cai(3)-Cai(1))./(2.*dx)).^2 + JCa(2)./Vnonjunct(2).*betai(2); 
dCaidt(3) = betai(3) .* (DCa + gammai(3).*DCaBm) .* ( (Cai(4)-2.*Cai(3)+Cai(2))./dx.^2 + (Cai(4)-Cai(2))./(2.*j(3).*dx.^2) ) - 2.*betai(3).*gammai(3).*DCaBm./(KdBCa + Cai(3)) .* ((Cai(4)-Cai(2))./(2.*dx)).^2 + JCa(3)./Vnonjunct(3).*betai(3); 
dCaidt(4) = betai(4) .* (DCa + gammai(4).*DCaBm) .* ( (Cai(4)-2.*Cai(4)+Cai(3))./dx.^2 + (Cai(4)-Cai(3))./(2.*j(4).*dx.^2) ) - 2.*betai(4).*gammai(4).*DCaBm./(KdBCa + Cai(4)) .* ((Cai(4)-Cai(3))./(2.*dx)).^2 + JCa(4)./Vnonjunct(4).*betai(4); 

dy(i_Cacenter:length(Cai)+i_Cacenter-1) = dCaidt;

dCaSRdt = zeros(length(CaSR),1); 
dCaSRdt(1) = betaSR(1) .* (DCaSR) .* ( (CaSR(2)-2.*CaSR(1)+CaSR(1))./dx.^2 + (CaSR(2)-CaSR(1))./(2.*j(1).*dx.^2) ) + JSRCa(1)./VSR(1).*betaSR(1); 
dCaSRdt(2) = betaSR(2) .* (DCaSR) .* ( (CaSR(3)-2.*CaSR(2)+CaSR(1))./dx.^2 + (CaSR(3)-CaSR(1))./(2.*j(2).*dx.^2) ) + JSRCa(2)./VSR(2).*betaSR(2); 
dCaSRdt(3) = betaSR(3) .* (DCaSR) .* ( (CaSR(4)-2.*CaSR(3)+CaSR(2))./dx.^2 + (CaSR(4)-CaSR(2))./(2.*j(3).*dx.^2) ) + JSRCa(3)./VSR(3).*betaSR(3); 
dCaSRdt(4) = betaSR(4) .* (DCaSR) .* ( (CaSR(4)-2.*CaSR(4)+CaSR(3))./dx.^2 + (CaSR(4)-CaSR(3))./(2.*j(4).*dx.^2) ) + JSRCa(4)./VSR(4).*betaSR(4); 

dy(i_Cacenter+length(j):length(j).*2+i_Cacenter-1) = dCaSRdt;

%% Output ion currents and fluxes
currents = [Iion Istim INa ICaL It Isus IKr IKs IK1 If INab ICab ICaP INaK INaCa i_KACh];
%currents = [Iion Istim INa ICaL It Isus IKr IKs IK1 If INab ICab ICaP INaK INaCa Jrelss Jrel1 Jrel2 Jrel3 J_bulkSERCA1 J_bulkSERCA2 J_bulkSERCA3 J_bulkSERCAss JSRCaleak1 JSRCaleak2 JSRCaleak3 JSRCaleakss];
%currents = [IKr IKs It i_KACh];

end % EOF