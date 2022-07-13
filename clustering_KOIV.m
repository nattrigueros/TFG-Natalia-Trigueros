%% CLUSTERING OF AP SIGNALS
% Multisignal 1-D Wavelet Analysis, a 1-D multisignal is a set of 1-D 
% signals of same length stored as a matrix organized rowwise (or 
% columnwise).

%% STEP 1: Pre-processing

% Load the appropiate stored .mat with the population in physiological
% ranges, for Koivumaki.
f=load (['PopPhysKOIV_control.mat']);
t=f.t;

g=load (['PopulationKOIV_raw.mat']);
s=g.s;

h=load (['NonPhysKOIV_AP.mat']);
nonphy=h.nonphy;

% Creation of a matrix X composed of:
% - Rows: Each member of the population
% - Columns: Each time instant
% The matrix is formed by the voltage value for each POM member at each 
% time instant.

% Creation of a matrix A composed of:
% - Rows: Each member of the population
% - Unique column: Each value for APD90
% The matrix is formed by the APD90 value for each POM member.

% Creation of X matrix
min=10000;
for i=1:length(t.result)
    l=length(t.result{i, 2}{1, 1}(:,1));
    if l<min
        min=l;
    end
end

x=zeros(length(t.result),min);

for i=1:length(t.result)
    h= t.result{i, 2}{1, 1}(1:min,1);
    h_2=transpose(h);
    x = [x;h_2];
end

% Creation of A matrix
a=zeros(length(t.result),1);

for i=1:length(t.result)
    h= t.biomarkers(i,3);
    a = [a;h];
end
       
X = x(length(t.result)+1:end,:); % Final signal matrix!
A = a(length(t.result)+1:end,:); % Final APD90 matrix!

[nbSIG,nbVAL] = size(X);

% Plot a graph with all the APs in the studied POM and the discarded ones
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for i=1:length(s.result)
    plot((s.result{i, 1}{1, 1}),(s.result{i, 2}{1, 1}(:,1)),'b')
    hold on
end

axis tight
axis([-0.01 0.4 -90 50])
xticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4])
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlabel(['Time (ms)'])
ylabel(['AP (mV)'])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')

% Plot a graph with all the APs in the studied POM
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for i=1:length(t.result)
    plot((t.result{i, 1}{1, 1}),(t.result{i, 2}{1, 1}(:,1)),'r')
    hold on
end

axis tight
axis([-0.01 0.4 -90 50])
xticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4])
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlabel(['Time (ms)'])
ylabel(['AP (mV)'])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')

% % Plot a graph with all the APs discarded from the POM
% figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
% 
% for i=1:length(nonphy.index)
%     ind=nonphy.index(i);
%     plot((s.result{ind, 1}{1, 1}),(s.result{ind, 2}{1, 1}(:,1)),'b')
%     hold on
% end
% 
% axis tight
% axis([-0.1 0.4 -90 40])
% xlabel(['Time (ms)'])
% ylabel(['AP (mV)'])

% % Plot a graph with a characteristic AP
% figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
% 
% plot((t.result{50, 1}{1, 1}),(t.result{50, 2}{1, 1}(:,1)),'b','LineWidth',1.7)
% hold on
% 
% axis tight
% axis([-0.01 0.4 -90 50])
% xticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4])
% xticklabels({'0','50','100','150','200','250','300','350','400'})
% xlabel(['Time (ms)'])
% ylabel(['AP (mV)'])


% Plot a graph with the choosen APs discarded from the POM
% Criteria followed for election was based on biomarkers values (the most
% noticeable).
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])


plot((t.result{50, 1}{1, 1}),(t.result{50, 2}{1, 1}(:,1)),'k','LineWidth',2.5)
hold on
plot((s.result{521, 1}{1, 1}),(s.result{521, 2}{1, 1}(:,1)),'m','LineWidth',1.7)
hold on
plot((s.result{693, 1}{1, 1}),(s.result{693, 2}{1, 1}(:,1)),'y','LineWidth',1.7)
hold on
plot((s.result{448, 1}{1, 1}),(s.result{448, 2}{1, 1}(:,1)),'c','LineWidth',1.7)
hold on

axis tight
axis([-0.01 0.4 -90 40])
xticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4])
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlabel(['Time (ms)'])
xlabel(['Time (s)'])
ylabel(['AP (mV)'])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')

%% STEP 2: Multisignal Row Decomposition
% Perform a wavelet decomposition at level 7 using the ' _sym4_ ' wavelet.
dirDec = 'r';         % Direction of decomposition
level  = 7;           % Level of decomposition  
wname  = 'sym4';      % Near symmetric wavelet
decROW = mdwtdec(dirDec,X,level,wname);

%% STEP 3: Clustering Row Signals
% Clustering offers a convenient procedure to summarize a large set of
% signals using sparse wavelet representations. We will implement
% hierarchical clustering using |mdwtcluster|

% Set to the number of clusters to 9 (nº of atrial regions).

P2 = mdwtcluster(decROW,'lst2clu',{'s'},'maxclust',9);

IdxInCluster = cell(1,9);
for k = 1:9
    IdxInCluster{k} = find(P2.IdxCLU==k);
end

% Perfom a plot of all the physiological APs in clusters with color distinctions.
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};

    if k==1 % BBLA
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 2}{1, 1}(:,1)),'color',[0.4940 0.1840 0.5560],'LineWidth',1.7)
        axis tight
        
    end
    if k==2 % LA
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 2}{1, 1}(:,1)),'y','LineWidth',1.7)
        axis tight
        
    end
    if k==3 % RA/PM
        hold on
        plot(t.result{idxK(10), 1}{1, 1},(t.result{idxK(10), 2}{1, 1}(:,1)),'b','LineWidth',1.7)
        axis tight
        
    end
    if k==4 % CT/BBRA
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 2}{1, 1}(:,1)),'c','LineWidth',1.7)
        axis tight
        
    end
    if k==5 % LAA
        hold on
        plot(t.result{idxK(3), 1}{1, 1},(t.result{idxK(3), 2}{1, 1}(:,1)),'color',[0.6350 0.0780 0.1840],'LineWidth',1.5)
        axis tight
        
    end
    if k==6 % PV
        hold on
        plot(t.result{idxK(19), 1}{1, 1},(t.result{idxK(19), 2}{1, 1}(:,1)),'g','LineWidth',1.5)
        axis tight
        
    end
    if k==7 % RAA
        hold on
        plot(t.result{idxK(4), 1}{1, 1},(t.result{idxK(4), 2}{1, 1}(:,1)),'color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
        axis tight
        
    end
    if k==8 % MVR
        hold on
        plot(t.result{idxK(16), 1}{1, 1},(t.result{idxK(16), 2}{1, 1}(:,1)),'m','LineWidth',1.5)
        axis tight
        
    end
    if k==9 % TVR
        hold on
        plot(t.result{idxK(end-1), 1}{1, 1},(t.result{idxK(end-1), 2}{1, 1}(:,1)),'r','LineWidth',1.5)
        axis tight
        axis([-0.01 0.4 -90 40])
        xticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4])
        xticklabels({'0','50','100','150','200','250','300','350','400'})
        xlabel(['Time (ms)'])
        ylabel(['AP (mV)'])
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
    end
    lgd=legend('LA','LAA','CT/BBRA','RAA','RA/PM','BBLA','MVR','TVR','PV');
    lgd.FontSize = 16;
    lgd.Title.String = 'Atrial Regions';
end 

%% STEP 4: Clustering Analysis
% - Find out average APD90 value for each cluster. (For posterior atrial
%   classification purposes).
% - Save the scaling factors set for each population member in the left
%   atrial regions for simulation in pAF.
% - Find out ONE representative AP for each cluster,
%   and perform plot.

%  REGION CLASSIFICATION based on APD90:
%  K = 1 --- BBLA
%  K = 2 --- LA
%  K = 3 --- RA/PM
%  K = 4 --- CT/BBRA
%  K = 5 --- LAA
%  K = 6 --- PV
%  K = 7 --- RAA
%  K = 8 --- MVR
%  K = 9 --- TVR

s_pv=zeros(25,16); % K=1 
s_bbla=zeros(23,16); % K=2 
s_laa=zeros(39,16); % K=6 
s_mvr=zeros(33,16); % K=8 
s_la=zeros(26,16); % K=9 

for k = 1:9 
    idxK = IdxInCluster{k};
    x=zeros(length(idxK),1);
    a=zeros(length(idxK),1); 

    if k==1 % LA
        hold on
        
        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k1 = mean(a);  

        % Save the scaling factors for BBLA region
        for j=1:length(idxK)
            s_la =[s_la;t.params{idxK(j),1}(1,:)];
        end
        s_la = s_la(length(idxK)+1:end,:);

    end

    if k==2 % LAA
        hold on

        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k2 = mean(a);

        % Save the scaling factors for LA region
        for j=1:length(idxK)
            s_laa =[s_laa;t.params{idxK(j),1}(1,:)];
        end
        s_laa = s_laa(length(idxK)+1:end,:);

    end

    if k==3  % CT/BBRA
        hold on

        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k3 = mean(a);
    end

    if k==4  % RAA
        hold on

        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k4 = mean(a);
    end

    if k==5  % RA/PM
        hold on

        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k5 = mean(a);

       
    end

    if k==6  % BBLA
        hold on

        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k6 = mean(a);

        % Save the scaling factors for LAA region
        for j=1:length(idxK)
            s_bbla =[s_bbla;t.params{idxK(j),1}(1,:)];
        end
        s_bbla = s_bbla(length(idxK)+1:end,:);

    end

    if k==7  % MVR
        hold on

        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k7 = mean(a);


        % Save the scaling factors for PV region
        for j=1:length(idxK)
            s_mvr =[s_mvr;t.params{idxK(j),1}(1,:)];
        end
        s_mvr = s_mvr(length(idxK)+1:end,:);
    end

    if k==8  % TVR
        hold on

        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k8 = mean(a);

    end

    if k==9  % PV
        hold on

        % Extract mean APD for this k POM members
        for j=1:length(idxK)
            a =[a;A(idxK(j))];
        end
        a = a(length(idxK)+1:end,:)
        mean_k9 = mean(a);

        % Save the scaling factors for PV region
        for j=1:length(idxK)
            s_pv =[s_pv;t.params{idxK(j),1}(1,:)];
        end
        s_pv = s_pv(length(idxK)+1:end,:);
    end
end

% Display table with mean APD90 for each k cluster
    disp(['APD90 mean values for each atrial cluster.'])
    VarNames = {'K=1', 'K=2', 'K=3', 'K=4','K=5', 'K=6', 'K=7', 'K=8', 'K=9'};
    apd_90_mean = table(mean_k1,mean_k2,mean_k3,mean_k4,mean_k5,mean_k6,mean_k7,mean_k8,mean_k9, 'VariableNames',VarNames)


% Perform a plot showing the APs that get in each cluster with different
% colors.
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    
    if k==6
        hold on
        subplot(3,3,2)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'g')
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 2}{1, 1}(:,1)),'k','LineWidth',1.5)
        
        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('BBLA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==1
        hold on
        subplot(3,3,9)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'color',[0.4940 0.1840 0.5560])
            hold on
        end
        hold on       
        plot(t.result{idxK(12), 1}{1, 1},(t.result{idxK(12), 2}{1, 1}(:,1)),'k','LineWidth',1.5)

        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('LA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
   
    if k==4
        hold on
        subplot(3,3,4)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'c')
            hold on
        end
        hold on
        plot(t.result{idxK(10), 1}{1, 1},(t.result{idxK(10), 2}{1, 1}(:,1)),'k','LineWidth',1.5)
        
        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('RAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

     if k==7
        hold on
        subplot(3,3,8)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'color',[0.9290 0.6940 0.1250])
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 2}{1, 1}(:,1)),'k','LineWidth',1.5)
        
        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('MVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
     end

     if k==9
        hold on
        subplot(3,3,1)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'r')
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 2}{1, 1}(:,1)),'k','LineWidth',1.5)
        
        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('PV')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
     end

     if k==5
        hold on
        subplot(3,3,7)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'color',[0.6350 0.0780 0.1840])
            hold on
        end
        hold on
        plot(t.result{idxK(7), 1}{1, 1},(t.result{idxK(7), 2}{1, 1}(:,1)),'k','LineWidth',1.5)
        
        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('RA/PM')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
     end

     if k==3
        hold on
        subplot(3,3,3)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'b')
            hold on
        end
        hold on
        plot(t.result{idxK(end-3), 1}{1, 1},(t.result{idxK(end-3), 2}{1, 1}(:,1)),'k','LineWidth',1.5)
        
        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('CT/BBRA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
     end

     if k==8
        hold on
        subplot(3,3,5)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'m')
            hold on
        end
        hold on
        plot(t.result{idxK(5), 1}{1, 1},(t.result{idxK(5), 2}{1, 1}(:,1)),'k','LineWidth',1.5)
        
        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('TVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==2
        hold on
        subplot(3,3,6)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 2}{1, 1}(:,1)),'y')
            hold on
        end
        hold on
        plot(t.result{idxK(end-8), 1}{1, 1},(t.result{idxK(end-8), 2}{1, 1}(:,1)),'k','LineWidth',1.5)
        
        axis([-0.01 0.4 -90 50])
        xlabel(['Time (ms)' '          (' int2str(length(idxK)) ')'])
        ylabel(['AP (mV)'])
        title('LAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -90 40])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
     
end 

%% STEP 5: Current study

% Sodium current

figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    
    if k==1 
        hold on
        subplot(3,3,9)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'color',[0.4940 0.1840 0.5560])
            hold on
        end
        hold on
        plot(t.result{idxK(12), 1}{1, 1},(t.result{idxK(12), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('LA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end

    if k==2 
        hold on
        subplot(3,3,6)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'y')
            hold on
        end
        hold on
        plot(t.result{idxK(end-8), 1}{1, 1},(t.result{idxK(end-8), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('LAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end

    if k==3 
        hold on
        subplot(3,3,3)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'b')
            hold on
        end
        hold on
        plot(t.result{idxK(end-3), 1}{1, 1},(t.result{idxK(end-3), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('CT/BBRA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end

    if k==4 
        hold on
        subplot(3,3,4)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'c')
            hold on
        end
        hold on
        plot(t.result{idxK(10), 1}{1, 1},(t.result{idxK(10), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('RAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end

    if k==5 
        hold on
        subplot(3,3,7)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'color',[0.6350 0.0780 0.1840])
            hold on
        end
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('RA/PM')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end
    if k==6 
        hold on
        subplot(3,3,2)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'g')
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('BBLA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end
    if k==7 
        hold on
        subplot(3,3,8)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'color',[0.9290 0.6940 0.1250])
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('MVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end
    if k==8 
        hold on
        subplot(3,3,5)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'m')
            hold on
        end
        hold on
        plot(t.result{idxK(28), 1}{1, 1},(t.result{idxK(28), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('TVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end
    if k==9 % 
        hold on
        subplot(3,3,1)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,3)),'r')
            hold on
        end
        
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,3)),'k','LineWidth',1.5)
        axis([0 0.01 -10000 30])  
        ylabel(['I (pA)'])
        title('PV')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        
        xticks([0 0.005 0.01])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end
    
end

     

% Acetylcholine current 

figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    
    if k==1 
        hold on
        subplot(3,3,9)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'color',[0.4940 0.1840 0.5560])
            hold on
        end
        hold on
        plot(t.result{idxK(12), 1}{1, 1},(t.result{idxK(12), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('LA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==2 
        hold on
        subplot(3,3,6)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'y')
            hold on
        end
        hold on
        plot(t.result{idxK(end-8), 1}{1, 1},(t.result{idxK(end-8), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('LAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==3 
        hold on
        subplot(3,3,3)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'b')
            hold on
        end
        hold on
        plot(t.result{idxK(end-3), 1}{1, 1},(t.result{idxK(end-3), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('CT/BBRA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==4 
        hold on
        subplot(3,3,4)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'c')
            hold on
        end
        hold on
        plot(t.result{idxK(10), 1}{1, 1},(t.result{idxK(10), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('RAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==5 
        hold on
        subplot(3,3,7)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'color',[0.6350 0.0780 0.1840])
            hold on
        end
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('RA/PM')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==6 
        hold on
        subplot(3,3,2)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'g')
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('BBLA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==7 
        hold on
        subplot(3,3,8)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'color',[0.9290 0.6940 0.1250])
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('MVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==8 
        hold on
        subplot(3,3,5)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'m')
            hold on
        end
        hold on
        plot(t.result{idxK(28), 1}{1, 1},(t.result{idxK(28), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('TVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==9 % 
        hold on
        subplot(3,3,1)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,16)),'r')
            hold on
        end
        
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,16)),'k','LineWidth',1.5)

        ylabel(['I (pA)'])
        title('PV')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -0 0.25])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
        xlabel(['Time (ms)'])
    end
    
end


% K1 current   

figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    
    if k==1 
        hold on
        subplot(3,3,9)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'color',[0.4940 0.1840 0.5560])
            hold on
        end
        hold on
        plot(t.result{idxK(12), 1}{1, 1},(t.result{idxK(12), 3}{1, 1}(:,9)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('LA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==2 
        hold on
        subplot(3,3,6)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'y')
            hold on
        end
        hold on
        plot(t.result{idxK(end-8), 1}{1, 1},(t.result{idxK(end-8), 3}{1, 1}(:,9)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('LAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==3 
        hold on
        subplot(3,3,3)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'b')
            hold on
        end
        hold on
        plot(t.result{idxK(end-3), 1}{1, 1},(t.result{idxK(end-3), 3}{1, 1}(:,9)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('CT/BBRA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==4 
        hold on
        subplot(3,3,4)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'c')
            hold on
        end
        hold on
        plot(t.result{idxK(10), 1}{1, 1},(t.result{idxK(10), 3}{1, 1}(:,9)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('RAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==5 
        hold on
        subplot(3,3,7)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'color',[0.6350 0.0780 0.1840])
            hold on
        end
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,9)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('RA/PM')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==6 
        hold on
        subplot(3,3,2)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'g')
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,9)),'k','LineWidth',1.5)
      
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('BBLA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==7 
        hold on
        subplot(3,3,8)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'color',[0.9290 0.6940 0.1250])
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,9)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('MVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==8 
        hold on
        subplot(3,3,5)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'m')
            hold on
        end
        hold on
        plot(t.result{idxK(28), 1}{1, 1},(t.result{idxK(28), 3}{1, 1}(:,9)),'k','LineWidth',1.5)

        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('TVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==9 % 
        hold on
        subplot(3,3,1)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,9)),'r')
            hold on
        end
        
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,9)),'k','LineWidth',1.5)

        ylabel(['I (pA)'])
        title('PV')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([-0.01 0.4 -10 60])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
        xlabel(['Time (ms)'])
    end
    
end

% Kr current  
            
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    
    if k==1 
        hold on
        subplot(3,3,9)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'color',[0.4940 0.1840 0.5560])
            hold on
        end
        hold on
        plot(t.result{idxK(12), 1}{1, 1},(t.result{idxK(12), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
       
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('LA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==2 
        hold on
        subplot(3,3,6)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'y')
            hold on
        end
        hold on
        plot(t.result{idxK(end-8), 1}{1, 1},(t.result{idxK(end-8), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
     
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('LAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==3 
        hold on
        subplot(3,3,3)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'b')
            hold on
        end
        hold on
        plot(t.result{idxK(end-3), 1}{1, 1},(t.result{idxK(end-3), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('CT/BBRA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==4 
        hold on
        subplot(3,3,4)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'c')
            hold on
        end
        hold on
        plot(t.result{idxK(10), 1}{1, 1},(t.result{idxK(10), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('RAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end

    if k==5 
        hold on
        subplot(3,3,7)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'color',[0.6350 0.0780 0.1840])
            hold on
        end
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('RA/PM')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==6 
        hold on
        subplot(3,3,2)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'g')
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('BBLA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==7 
        hold on
        subplot(3,3,8)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'color',[0.9290 0.6940 0.1250])
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('MVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==8 
        hold on
        subplot(3,3,5)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'m')
            hold on
        end
        hold on
        plot(t.result{idxK(28), 1}{1, 1},(t.result{idxK(28), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('TVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
    end
    if k==9 % 
        hold on
        subplot(3,3,1)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,7)),'r')
            hold on
        end
        
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,7)),'k','LineWidth',1.5)
        
        ylabel(['I (pA)'])
        title('PV')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.4 0 4.5])
        xticks([0 0.2 0.4])
        xticklabels({'0','200','400'})
        xlabel(['Time (ms)'])
    end
    
end

% CaL current

figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    
    if k==1 
        hold on
        subplot(3,3,9)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'color',[0.4940 0.1840 0.5560])
            hold on
        end
        hold on
        plot(t.result{idxK(12), 1}{1, 1},(t.result{idxK(12), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('LA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
    end

    if k==2 
        hold on
        subplot(3,3,6)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'y')
            hold on
        end
        hold on
        plot(t.result{idxK(end-8), 1}{1, 1},(t.result{idxK(end-8), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('LAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
    end

    if k==3 
        hold on
        subplot(3,3,3)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'b')
            hold on
        end
        hold on
        plot(t.result{idxK(end-3), 1}{1, 1},(t.result{idxK(end-3), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('CT/BBRA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
    end

    if k==4 
        hold on
        subplot(3,3,4)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'c')
            hold on
        end
        hold on
        plot(t.result{idxK(10), 1}{1, 1},(t.result{idxK(10), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('RAA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
    end

    if k==5 
        hold on
        subplot(3,3,7)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'color',[0.6350 0.0780 0.1840])
            hold on
        end
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('RA/PM')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
    end
    if k==6 
        hold on
        subplot(3,3,2)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'g')
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('BBLA')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
    end
    if k==7 
        hold on
        subplot(3,3,8)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'color',[0.9290 0.6940 0.1250])
            hold on
        end
        hold on
        plot(t.result{idxK(1), 1}{1, 1},(t.result{idxK(1), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('MVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
    end
    if k==8 
        hold on
        subplot(3,3,5)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'m')
            hold on
        end
        hold on
        plot(t.result{idxK(28), 1}{1, 1},(t.result{idxK(28), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        xlabel(['Time (ms)'])
        ylabel(['I (pA)'])
        title('TVR')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
    end
    if k==9 % 
        hold on
        subplot(3,3,1)
        for i=1:length(idxK)
            plot((t.result{idxK(i), 1}{1, 1}),(t.result{idxK(i), 3}{1, 1}(:,4)),'r')
            hold on
        end
        
        hold on
        plot(t.result{idxK(8), 1}{1, 1},(t.result{idxK(8), 3}{1, 1}(:,4)),'k','LineWidth',1.5)
        
        ylabel(['I (pA)'])
        title('PV')
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
        set(gca,'XTickLabelMode','auto')
        set(gca,'XTick',[])
        axis([0 0.1 -600 5])
        xticks([0 0.05 0.1])
        xticklabels({'0','50','100'})
        xlabel(['Time (ms)'])
    end
    
end
%% STEP 6: Save scaling factors for Left Atrial Regions.

pv.name='ScalingFactors_PV';
pv.params=s_pv;
save(['ScalingFactorsKOIV_PV.mat'], 'pv')

bbla.name='ScalingFactors_BBLA';
bbla.params=s_bbla;
save(['ScalingFactorsKOIV_BBLA.mat'], 'bbla')

laa.name='ScalingFactors_LAA';
laa.params=s_laa;
save(['ScalingFactorsKOIV_LAA.mat'], 'laa')

mvr.name='ScalingFactors_MVR';
mvr.params=s_mvr;
save(['ScalingFactorsKOIV_MVR.mat'], 'mvr')

la.name='ScalingFactors_LA';
la.params=s_la;
save(['ScalingFactorsKOIV_LA.mat'], 'la')

%% STEP 7: See conductivity trends in each Atrial Region

% Gto
gto_data=zeros(289,1);
k_data=zeros(289,1);
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    if k==1
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end 
    if k==2
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==3
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==4
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==5
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==6
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==7
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==8
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==9
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(5);
            gto_data = [gto_data;data];
            k_data = [k_data;k];
        end 
    end
end 
Gto_data = gto_data(289+1:end,:);
K_data = k_data(289+1:end,:);
gscatter(1:289,Gto_data,K_data)
hold on
plot(1:289,ones([1,289]),'k','LineWidth',2)
lgd=legend('LA','LAA','CT/BBRA','RAA','RA/PM','BBLA','MVR','TVR','PV');
lgd.FontSize = 14;
lgd.Location='northeast';
lgd.NumColumns=3;
lgd.Title.String = 'Atrial Regions';
ylabel(['Multiplication factor'])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')
set(gca,'XTick',[])

% Gcal
gcal_data=zeros(289,1);
k_data=zeros(289,1);
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    if k==1
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end 
    if k==2
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==3
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==4
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==5
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==6
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==7
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==8
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==9
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(4);
            gcal_data = [gcal_data;data];
            k_data = [k_data;k];
        end 
    end
end 
Gcal_data = gcal_data(289+1:end,:);
K_data = k_data(289+1:end,:);
gscatter(1:289,Gcal_data,K_data)
hold on
plot(1:289,ones([1,289]),'k','LineWidth',2)
lgd=legend('LA','LAA','CT/BBRA','RAA','RA/PM','BBLA','MVR','TVR','PV');
lgd.FontSize = 14;
lgd.Location='northeast';
lgd.NumColumns=3;
lgd.Title.String = 'Atrial Regions';
ylabel(['Multiplication factor'])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')
set(gca,'XTick',[])

% Gkr
gkr_data=zeros(289,1);
k_data=zeros(289,1);
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    if k==1
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end 
    if k==2
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==3
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==4
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==5
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==6
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==7
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==8
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==9
        for i=1:length(idxK)
            data= t.params{idxK(i), 1}(7);
            gkr_data = [gkr_data;data];
            k_data = [k_data;k];
        end 
    end
end 
Gkr_data = gkr_data(289+1:end,:);
K_data = k_data(289+1:end,:);
gscatter(1:289,Gkr_data,K_data)
hold on
plot(1:289,ones([1,289]),'k','LineWidth',2)
lgd=legend('LA','LAA','CT/BBRA','RAA','RA/PM','BBLA','MVR','TVR','PV');
lgd.FontSize = 14;
lgd.Location='northeast';
lgd.NumColumns=3;
lgd.Title.String = 'Atrial Regions';
ylabel(['Multiplication factor'])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')
set(gca,'XTick',[])

% APD
apd_data=zeros(289,1);
k_data=zeros(289,1);
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for k = 1:9
    idxK = IdxInCluster{k};
    if k==1
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end 
    if k==2
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==3
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==4
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==5
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==6
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==7
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==8
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end
    if k==9
        for i=1:length(idxK)
            data= t.biomarkers(idxK(i),3);
            apd_data = [apd_data;data];
            k_data = [k_data;k];
        end 
    end
    
end 
APD_data = apd_data(289+1:end,:);
K_data = k_data(289+1:end,:);
gscatter(1:289,APD_data,K_data)
hold on
plot(1:289,ones([1,289])*233.02,'k','LineWidth',2)
lgd=legend('LA','LAA','CT/BBRA','RAA','RA/PM','BBLA','MVR','TVR','PV');
lgd.FontSize = 14;
lgd.Location='northeast';
lgd.NumColumns=3;
lgd.Title.String = 'Atrial Regions';
ylabel(['APD 90 (ms)'])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')
set(gca,'XTick',[])