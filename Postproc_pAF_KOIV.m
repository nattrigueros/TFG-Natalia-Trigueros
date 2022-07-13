%% POST-PROCESSING WITH pAF REMODELLING

%% STEP 1: Pre-processing

% Load the appropiate stored .mat with the population in physiological
% ranges, for the five LA regions.

f=load (['PopPhysKOIV_pAF_PV.mat']);
PV=f.t;

f=load (['PopPhysKOIV_pAF_BBLA.mat']);
BBLA=f.t;

f=load (['PopPhysKOIV_pAF_LAA.mat']);
LAA=f.t;

f=load (['PopPhysKOIV_pAF_MVR.mat']);
MVR=f.t;

f=load (['PopPhysKOIV_pAF_LA.mat']);
LA=f.t;

% Plot a graph with all the APs in the studied POM
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for i=1:length(PV.result)
    plot((PV.result{i, 1}{1, 1}),(PV.result{i, 2}{1, 1}(:,1)),'r')
    hold on
end
for i=1:length(BBLA.result)
    plot((BBLA.result{i, 1}{1, 1}),(BBLA.result{i, 2}{1, 1}(:,1)),'r')
    hold on
end
for i=1:length(LAA.result)
    plot((LAA.result{i, 1}{1, 1}),(LAA.result{i, 2}{1, 1}(:,1)),'r')
    hold on
end
for i=1:length(MVR.result)
    plot((MVR.result{i, 1}{1, 1}),(MVR.result{i, 2}{1, 1}(:,1)),'r')
    hold on
end
for i=1:length(LA.result)
    plot((LA.result{i, 1}{1, 1}),(LA.result{i, 2}{1, 1}(:,1)),'r')
    hold on
end
axis tight
axis([-0.01 0.4 -90 50])

xlabel(['Time (ms)'])
ylabel(['AP (mV)'])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')
xticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4])
xticklabels({'0','50','100','150','200','250','300','350','400'})

% Plot with just representative APs
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
plot((PV.result{5, 1}{1, 1}),(PV.result{5, 2}{1, 1}(:,1)),'r','LineWidth',1.7)
hold on

plot((BBLA.result{4, 1}{1, 1}),(BBLA.result{4, 2}{1, 1}(:,1)),'g','LineWidth',1.7)
hold on

plot((LAA.result{8, 1}{1, 1}),(LAA.result{8, 2}{1, 1}(:,1)),'y','LineWidth',1.7)
hold on


plot((MVR.result{16, 1}{1, 1}),(MVR.result{16, 2}{1, 1}(:,1)),'color',[0.9290 0.6940 0.1250],'LineWidth',1.7)
hold on

plot((LA.result{6, 1}{1, 1}),(LA.result{6, 2}{1, 1}(:,1)),'color',[0.4940 0.1840 0.5560],'LineWidth',1.7)
hold on

axis tight
axis([-0.01 0.4 -90 50])

xlabel(['Time (ms)'])
ylabel(['AP (mV)'])
lgd=legend('PV','BBLA','LAA','MVR','LA')
lgd.FontSize = 16;
lgd.Title.String = 'Atrial Regions';
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
set(gca,'XTickLabelMode','auto')
xticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4])
xticklabels({'0','50','100','150','200','250','300','350','400'})

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

regions = [PV,BBLA,LAA,MVR,LA];
s_pv=zeros(6,16); 
s_bbla=zeros(11,16); 
s_laa=zeros(33,16); 
s_mvr=zeros(25,16);
s_la=zeros(19,16);

s_pv_short=zeros(1,16); 
s_bbla_short=zeros(1,16); 
s_laa_short=zeros(1,16); 
s_mvr_short=zeros(1,16);
s_la_short=zeros(1,16);

s_pv_long=zeros(1,16); 
s_bbla_long=zeros(1,16); 
s_laa_long=zeros(1,16); 
s_mvr_long=zeros(1,16);
s_la_long=zeros(1,16);


for r=1:length(regions)
    min=10000;
    for i=1:length(regions(r).result)
        l=length(regions(r).result{i, 2}{1, 1}(:,1));
        if l<min
            min=l;
        end
    end
    
    x=zeros(length(regions(r).result),min);
    
    for i=1:length(regions(r).result)
        h= regions(r).result{i, 2}{1, 1}(1:min,1);
        h_2=transpose(h);
        x = [x;h_2];
    end
    
    % Creation of A matrix
    a=zeros(length(regions(r).result),1);
    
    for i=1:length(regions(r).result)
        h= regions(r).biomarkers(i,3);
        a = [a;h];
    end
           
    X = x(length(regions(r).result)+1:end,:); % Final signal matrix!
    A = a(length(regions(r).result)+1:end,:); % Final APD90 matrix!

    x=zeros(length(regions(r).result),1);
    a=zeros(length(regions(r).result),1); 

    if r==1  % PV
        hold on
        
        % Extract mean APD for this k POM members
        min1=10000;
        min2=0;
        for j=1:length(regions(r).result)
            a =[a;A(j)];
            if A(j)<min1
                min1=A(j);
                memb1=j;
            end
            if A(j)>min2
                min2=A(j);
                memb2=j;
            end
        end
        disp([ 'Member with lowest APD for PV is: ' int2str(memb1) ' with APD ' int2str(min1) ])
        disp(['Member with highest APD for PV is: ' int2str(memb2) ' with APD ' int2str(min2) ])

        a = a(length(regions(r).result)+1:end,:)
        mean_r1 = mean(a);  

        % Save the scaling factors for PV region
        for j=1:length(regions(r).result)
            s_pv =[s_pv;regions(r).params{j,1}(1,:)];
        end
        s_pv = s_pv(length(regions(r).result)+1:end,:);
        s_pv_short=regions(r).params{memb1,1}(1,:);
        s_pv_long=regions(r).params{memb2,1}(1,:);
    end

    if r==2  % BBLA
        hold on
        
        % Extract mean APD for this k POM members
        min1=10000;
        min2=0;
        for j=1:length(regions(r).result)
            a =[a;A(j)];
            if A(j)<min1
                min1=A(j);
                memb1=j;
            end
            if A(j)>min2
                min2=A(j);
                memb2=j;
            end
        end
        disp([ 'Member with lowest APD for BBLA is: ' int2str(memb1) ' with APD ' int2str(min1) ])
        disp(['Member with highest APD for BBLA is: ' int2str(memb2) ' with APD ' int2str(min2) ])

        a = a(length(regions(r).result)+1:end,:)
        mean_r2 = mean(a);  

        % Save the scaling factors for BBLA region
        for j=1:length(regions(r).result)
            s_bbla =[s_bbla;regions(r).params{j,1}(1,:)];
        end
        s_bbla = s_bbla(length(regions(r).result)+1:end,:);
        s_bbla_short=regions(r).params{memb1,1}(1,:);
        s_bbla_long=regions(r).params{memb2,1}(1,:);
    end

    if r==3  % LAA
        hold on
        
        % Extract mean APD for this k POM members
        min1=10000;
        min2=0;
        for j=1:length(regions(r).result)
            a =[a;A(j)];
            if A(j)<min1
                min1=A(j);
                memb1=j;
            end
            if A(j)>min2
                min2=A(j);
                memb2=j;
            end
        end
        disp([ 'Member with lowest APD for LAA is: ' int2str(memb1) ' with APD ' int2str(min1) ])
        disp(['Member with highest APD for LAA is: ' int2str(memb2) ' with APD ' int2str(min2) ])

        a = a(length(regions(r).result)+1:end,:)
        mean_r3 = mean(a);  

        % Save the scaling factors for LAA region
        for j=1:length(regions(r).result)
            s_laa =[s_laa;regions(r).params{j,1}(1,:)];
        end
        s_laa = s_laa(length(regions(r).result)+1:end,:);
        s_laa_short=regions(r).params{memb1,1}(1,:);
        s_laa_long=regions(r).params{memb2,1}(1,:);
    end
    
    if r==4  % MVR
        hold on
        
        % Extract mean APD for this k POM members
        min1=10000;
        min2=0;
        for j=1:length(regions(r).result)
            a =[a;A(j)];
            if A(j)<min1
                min1=A(j);
                memb1=j;
            end
            if A(j)>min2
                min2=A(j);
                memb2=j;
            end
        end
        disp([ 'Member with lowest APD for MVR is: ' int2str(memb1) ' with APD ' int2str(min1) ])
        disp(['Member with highest APD for MVR is: ' int2str(memb2) ' with APD ' int2str(min2) ])

        a = a(length(regions(r).result)+1:end,:)
        mean_r4 = mean(a);  

        % Save the scaling factors for PV region
        for j=1:length(regions(r).result)
            s_mvr =[s_mvr;regions(r).params{j,1}(1,:)];
        end
        s_mvr = s_mvr(length(regions(r).result)+1:end,:);
        s_mvr_short=regions(r).params{memb1,1}(1,:);
        s_mvr_long=regions(r).params{memb2,1}(1,:);
    end

    if r==5  % LA
        hold on
        
        % Extract mean APD for this k POM members
        min1=10000;
        min2=0;
        for j=1:length(regions(r).result)
            a =[a;A(j)];
            if A(j)<min1
                min1=A(j);
                memb1=j;
            end
            if A(j)>min2
                min2=A(j);
                memb2=j;
            end
        end
        disp([ 'Member with lowest APD for LA is: ' int2str(memb1) ' with APD ' int2str(min1) ])
        disp(['Member with highest APD for LA is: ' int2str(memb2) ' with APD ' int2str(min2) ])

        a = a(length(regions(r).result)+1:end,:)
        mean_r5 = mean(a);  

        % Save the scaling factors for PV region
        for j=1:length(regions(r).result)
            s_la =[s_la;regions(r).params{j,1}(1,:)];
        end
        s_la = s_la(length(regions(r).result)+1:end,:);
        s_la_short=regions(r).params{memb1,1}(1,:);
        s_la_long=regions(r).params{memb2,1}(1,:);
    end

end

% Display table with mean APD90 for each region
    disp(['APD90 mean values for each atrial cluster.'])
    VarNames = {'PV', 'BBLA', 'LAA', 'MVR','LA'};
    apd_90_mean = table(mean_r1,mean_r2,mean_r3,mean_r4,mean_r5, 'VariableNames',VarNames)

%% STEP 2: Current study.

% Sodium
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for r=1:length(regions)
for i=1:length(regions(r).result)
        if r==1 % PV
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,3)),'r')
            hold on
            plot(regions(r).result{5, 1}{1, 1},(regions(r).result{5, 3}{1, 1}(:,3)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('PV')
            axis([0 0.01 -10000 30])
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            
            xticks([0 0.005 0.01])
            xticklabels({'0','50','100'})
            xlabel(['Time (ms)'])
              
       
        end
        if r==2
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,3)),'g')
            hold on
            plot(regions(r).result{4, 1}{1, 1},(regions(r).result{4, 3}{1, 1}(:,3)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('BBLA')
            axis([0 0.01 -10000 30])
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            
            xticks([0 0.005 0.01])
            xticklabels({'0','50','100'})
            xlabel(['Time (ms)'])
        end
        if r==3
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,3)),'y')
            hold on
            plot(regions(r).result{8, 1}{1, 1},(regions(r).result{8, 3}{1, 1}(:,3)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('LAA')
            axis([0 0.01 -10000 30])
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            
            xticks([0 0.005 0.01])
            xticklabels({'0','50','100'})
            xlabel(['Time (ms)'])
        end
        if r==4
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,3)),'color',[0.9290 0.6940 0.1250])
            hold on
            plot(regions(r).result{16, 1}{1, 1},(regions(r).result{16, 3}{1, 1}(:,3)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('MVR')
            axis([0 0.01 -10000 30])
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            
            xticks([0 0.005 0.01])
            xticklabels({'0','50','100'})
            xlabel(['Time (ms)'])
        end
        if r==5
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,3)),'color',[0.4940 0.1840 0.5560])
            hold on
            plot(regions(r).result{6, 1}{1, 1},(regions(r).result{6, 3}{1, 1}(:,3)),'k','LineWidth',1.5)
            axis tight
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
end
end

% Acetylcholine
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for r=1:length(regions)
for i=1:length(regions(r).result)
        if r==1 % PV
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,16)),'r')
            hold on
            plot(regions(r).result{5, 1}{1, 1},(regions(r).result{5, 3}{1, 1}(:,16)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('PV')
            xlabel(['Time (ms)'])
            
            
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -0 0.35])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
        end
        if r==2
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,16)),'g')
            hold on
            plot(regions(r).result{4, 1}{1, 1},(regions(r).result{4, 3}{1, 1}(:,16)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('BBLA')
            
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -0 0.35])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==3
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,16)),'y')
            hold on
            plot(regions(r).result{8, 1}{1, 1},(regions(r).result{8, 3}{1, 1}(:,16)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('LAA')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -0 0.35])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==4
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,16)),'color',[0.9290 0.6940 0.1250])
            hold on
            plot(regions(r).result{16, 1}{1, 1},(regions(r).result{16, 3}{1, 1}(:,16)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('MVR')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -0 0.35])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==5
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,16)),'color',[0.4940 0.1840 0.5560])
            hold on
            plot(regions(r).result{6, 1}{1, 1},(regions(r).result{6, 3}{1, 1}(:,16)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('LA')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -0 0.35])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
end
end

% K1
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for r=1:length(regions)
for i=1:length(regions(r).result)
        if r==1 % PV
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,9)),'r')
            hold on
            plot(regions(r).result{5, 1}{1, 1},(regions(r).result{5, 3}{1, 1}(:,9)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('PV')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -10 70])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==2
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,9)),'g')
            hold on
            plot(regions(r).result{4, 1}{1, 1},(regions(r).result{4, 3}{1, 1}(:,9)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('BBLA')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -10 70])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==3
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,9)),'y')
            hold on
            plot(regions(r).result{8, 1}{1, 1},(regions(r).result{8, 3}{1, 1}(:,9)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('LAA')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -10 70])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==4
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,9)),'color',[0.9290 0.6940 0.1250])
            hold on
            plot(regions(r).result{16, 1}{1, 1},(regions(r).result{16, 3}{1, 1}(:,9)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('MVR')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -10 70])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==5
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,9)),'color',[0.4940 0.1840 0.5560])
            hold on
            plot(regions(r).result{6, 1}{1, 1},(regions(r).result{6, 3}{1, 1}(:,9)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('LA')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([-0.01 0.4 -10 70])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
end
end

% Kr
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for r=1:length(regions)
for i=1:length(regions(r).result)
        if r==1 % PV
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,7)),'r')
            hold on
            plot(regions(r).result{5, 1}{1, 1},(regions(r).result{5, 3}{1, 1}(:,7)),'k','LineWidth',1.5)
            axis tight
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
        if r==2
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,7)),'g')
            hold on
            plot(regions(r).result{4, 1}{1, 1},(regions(r).result{4, 3}{1, 1}(:,7)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('BBLA')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([0 0.4 0 4.5])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==3
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,7)),'y')
            hold on
            plot(regions(r).result{8, 1}{1, 1},(regions(r).result{8, 3}{1, 1}(:,7)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('LAA')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([0 0.4 0 4.5])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==4
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,7)),'color',[0.9290 0.6940 0.1250])
            hold on
            plot(regions(r).result{16, 1}{1, 1},(regions(r).result{16, 3}{1, 1}(:,7)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('MVR')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([0 0.4 0 4.5])
            xticks([0 0.2 0.4])
            xticklabels({'0','200','400'})
            xlabel(['Time (ms)'])
        end
        if r==5
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,7)),'color',[0.4940 0.1840 0.5560])
            hold on
            plot(regions(r).result{6, 1}{1, 1},(regions(r).result{6, 3}{1, 1}(:,7)),'k','LineWidth',1.5)
            axis tight
            ylabel(['I (pA)'])
            title('LA')
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
end

% CaL
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for r=1:length(regions)
for i=1:length(regions(r).result)
        if r==1 % PV
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,4)),'r')
            hold on
            plot(regions(r).result{5, 1}{1, 1},(regions(r).result{5, 3}{1, 1}(:,4)),'k','LineWidth',1.5)
            axis tight
            
            xlabel(['Time (ms)'])
            ylabel(['I (pA)'])
            title('PV')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',16) %'FontWeight','bold'
            set(gca,'XTickLabelMode','auto')
            set(gca,'XTick',[])
            axis([0 0.1 -600 5])
            xticks([0 0.05 0.1])
            xticklabels({'0','50','100'})
        end
        if r==2
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,4)),'g')
            hold on
            plot(regions(r).result{4, 1}{1, 1},(regions(r).result{4, 3}{1, 1}(:,4)),'k','LineWidth',1.5)
            axis tight
            
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
        if r==3
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,4)),'y')
            hold on
            plot(regions(r).result{8, 1}{1, 1},(regions(r).result{8, 3}{1, 1}(:,4)),'k','LineWidth',1.5)
            axis tight
            
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
        if r==4
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,4)),'color',[0.9290 0.6940 0.1250])
            hold on
            plot(regions(r).result{16, 1}{1, 1},(regions(r).result{16, 3}{1, 1}(:,4)),'k','LineWidth',1.5)
            axis tight
            
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
        if r==5
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 3}{1, 1}(:,4)),'color',[0.4940 0.1840 0.5560])
            hold on
            plot(regions(r).result{6, 1}{1, 1},(regions(r).result{6, 3}{1, 1}(:,4)),'k','LineWidth',1.5)
            axis tight
            
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
end
end
   
%% STEP 3: AP study.
figure('Units','normalized','Position',[0.2 0.2 0.6 0.6])
for r=1:length(regions)
     % Plot a graph with all the APs in the studied POM
    for i=1:length(regions(r).result)
        if r==1 % PV
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 2}{1, 1}(:,1)),'r')
            hold on
            plot(regions(r).result{5, 1}{1, 1},(regions(r).result{5, 2}{1, 1}(:,1)),'k','LineWidth',1.5)
            axis tight
            xlabel(['Time (ms)' '          (' int2str(length(regions(r).result)) ')'])
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
        if r==2
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 2}{1, 1}(:,1)),'g')
            hold on
            plot(regions(r).result{4, 1}{1, 1},(regions(r).result{4, 2}{1, 1}(:,1)),'k','LineWidth',1.5)
            axis tight
            xlabel(['Time (ms)' '          (' int2str(length(regions(r).result)) ')'])
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
        if r==3
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 2}{1, 1}(:,1)),'y')
            hold on
            plot(regions(r).result{8, 1}{1, 1},(regions(r).result{8, 2}{1, 1}(:,1)),'k','LineWidth',1.5)
            axis tight
            xlabel(['Time (ms)' '          (' int2str(length(regions(r).result)) ')'])
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
        if r==4
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 2}{1, 1}(:,1)),'color',[0.9290 0.6940 0.1250])
            hold on
            plot(regions(r).result{16, 1}{1, 1},(regions(r).result{16, 2}{1, 1}(:,1)),'k','LineWidth',1.5)
            axis tight
            xlabel(['Time (ms)' '          (' int2str(length(regions(r).result)) ')'])
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
        if r==5
            subplot(2,3,r)
            hold on
            plot((regions(r).result{i, 1}{1, 1}),(regions(r).result{i, 2}{1, 1}(:,1)),'color',[0.4940 0.1840 0.5560])
            hold on
            plot(regions(r).result{6, 1}{1, 1},(regions(r).result{6, 2}{1, 1}(:,1)),'k','LineWidth',1.5)
            axis tight
            xlabel(['Time (ms)' '          (' int2str(length(regions(r).result)) ')'])
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
    end
end 
    

%% STEP 4: Save scaling factors for longest and shortest APD90 in each Left Atrial Region.

pv.name='SF-APD_KOIV_PV';
pv.short=s_pv_short;
pv.long=s_pv_long;
save(['SF-APD_KOIV_PV.mat'], 'pv')

bbla.name='SF-APD_KOIV_BBLA';
bbla.short=s_bbla_short;
bbla.long=s_bbla_long;
save(['SF-APD_KOIV_BBLA.mat'], 'bbla')

laa.name='SF-APD_KOIV_LAA';
laa.short=s_laa_short;
laa.long=s_laa_long;
save(['SF-APD_KOIV_LAA.mat'], 'laa')

mvr.name='SF-APD_KOIV_MVR';
mvr.short=s_mvr_short;
mvr.long=s_mvr_long;
save(['SF-APD_KOIV_MVR.mat'], 'mvr')

la.name='SF-APD_KOIV_LA';
la.short=s_la_short;
la.long=s_la_long;
save(['SF-APD_KOIV_LA.mat'], 'la')

la_t.name='SF-t-APD_KOIV_LA';
la_t.sf=s_la;
save(['SF-t-APD_KOIV_LA.mat'], 'la_t')