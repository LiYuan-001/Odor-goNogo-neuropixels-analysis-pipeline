function [ISI,Explained,PC] = gngPCA_T1_T5_250411(odorIndex_allCells, rst_all_allCells, Time, x_range, numPC,meanFivePercentile,inhStimTag,mode,BehavPerformance, ISI_ABsession, xRange, yRange)

%calculates PCA for T1 - T5
% 200820 Kei

tic

bin=0.05;%sec, bin for psth
bin_size=50;%ms, bin for psth
win=[-2 6];%sec, bin for raster and psth(8 sec, 161 bins data)
showFig =1;
num_xrange=size(x_range, 2);
x_range2 = 1:1:num_xrange;
titleText = {'T1', 'T2', 'T3', 'T4', 'T5'};



fig=figure(1);
set( fig,'position',[20,20,1600,900]);
set(gcf,'renderer','Painters')


 

for ss= 1:5 % T1 - T5
    % calculate PSTH for each cell
    for i= 1:size(odorIndex_allCells,2) % cell 1 - cell 239
        
        odorIndex = odorIndex_allCells{1,i};
        
        
        rst_all = rst_all_allCells{1,i};
        for s= 1:4 % A, B, C, D
            trialID=find(odorIndex ==s);
            numTrials=size(trialID,1);
            increment =round((numTrials-10)/4);
            
            rst=[];
   
                
                trialRange(1,:) = 1:10; % first10
                
                %trialRange(2,:) = increment:increment+9; %T2
                    startT2 = round((5+numTrials/2)/2)-4;
                    endT2  =  round((5+numTrials/2)/2)+5;
                trialRange(2,:) = startT2:endT2  ; %T2
                
                trialRange(3,:) = round(numTrials/2)-4:round(numTrials/2)+5;  % middle10
                
                %trialRange(4,:) = increment*3:increment*3+9;  %T4
                    startT4 = round((numTrials/2 + numTrials-5)/2)-4;
                    endT4  =  round((numTrials/2 + numTrials-5)/2)+5;
                trialRange(4,:) = startT4:endT4  ; %T2    
                
                trialRange(5,:) = numTrials-9:numTrials; % last10
 

            
            for k=trialRange(ss,:)
                
                rst_oneTrial = rst_all(rst_all(:,2)==trialID(k),:) ;
                
                rst =[rst;rst_oneTrial];
            end
            
            psth = ana_peth_ss(rst,win,bin);
            psth_gs(i,:,s)=smoothing_gaussian_ss(psth,100,bin_size,2);
            
        end
    end
    
    tmpd1 = psth_gs(:,x_range,1);
    tmpd2 = psth_gs(:,x_range,2);
    tmpd3 = psth_gs(:,x_range,3);
    tmpd4 = psth_gs(:,x_range,4);
    CTS_data2=[tmpd1 tmpd2 tmpd3 tmpd4];
    
    CTS_data2=CTS_data2-min(CTS_data2,[],2);
    CTS_data2=CTS_data2./max(CTS_data2,[],2);
    CTS_data2=CTS_data2(~isnan(CTS_data2(:,1)),:); %removes cells with NaN data
    
    sigma = 200; %ms
    binSize = 50; %ms
    
    CTS_data2(:,x_range2) =   smooth_gaussian2(CTS_data2(:,x_range2),  sigma,binSize,1);
    CTS_data2(:,x_range2+num_xrange) = smooth_gaussian2(CTS_data2(:,x_range2+num_xrange),sigma,binSize,1);
    CTS_data2(:,x_range2+num_xrange*2) = smooth_gaussian2(CTS_data2(:,x_range2+num_xrange*2),sigma,binSize,1);
    CTS_data2(:,x_range2+num_xrange*3) = smooth_gaussian2(CTS_data2(:,x_range2+num_xrange*3),sigma,binSize,1);
    
    
    
    %% PCA
    [~,score,~,~, explained] = pca(CTS_data2');
    
    
    PC_A1=score(x_range2,1:numPC); %numPC =2 or 3
    PC_B1=score(x_range2+num_xrange,1:numPC);
    PC_C1=score(x_range2+num_xrange*2,1:numPC);
    PC_D1=score(x_range2+num_xrange*3,1:numPC);
    
%     PC.A1(ss) = PC_A1; %commented out on 7/10/23, assignment errors
%     PC.B1(ss) = PC_B1;
%     PC.C1(ss) = PC_C1;
%     PC.D1(ss) = PC_D1;
    
    Explained.PC1(ss)=sum(explained(1));
    Explained.PC2(ss)=sum(explained(2));
    if numPC>2
    Explained.PC3(ss)=sum(explained(3));
    end
    % ExplainedByPC1_PC6=sum(explained(1:6))
    % ExplainedByPC1_PC10=sum(explained(1:10))
    % numCells = numel(explained)
    
    
    % %%% Euclidean distance
    normalization = 0;
    % --------------------------
    ECdist = gngSessionSum2_PCA_ECdist_v5(PC_A1, PC_B1, PC_C1, PC_D1, normalization,4);
    
    % --------------------------
    
    
    %% calculate shuffling data
    
        % mean of 0.5-1.5 s after odor onset
        % 1-20 -1-0s 
        % 21-40 0-1s
        % 41-60 1-2s
        % 61-80 2-3s 
        % 81-100 3-4s 
        % 101-120 4-5s
    
    odor1000.AB(ss) =  mean(ECdist.AB1(31:50));
    odor1000.AC(ss) =  mean(ECdist.AC1(31:50));
    odor1000.AD(ss) =  mean(ECdist.AD1(31:50));
    odor1000.BC(ss) =  mean(ECdist.BC1(31:50));
    odor1000.BD(ss) =  mean(ECdist.BD1(31:50));
    odor1000.CD(ss) =  mean(ECdist.CD1(31:50));
    
    
    
    %% PCA continuous plot
    
    
    
    mMin=min(min(score));
    mMax=max(max(score));
    
    tp1=21;%0s 
    tp2=41;%1s
    tp3=61;%2s
    tp4=81;%3s
    tp5=101;%4s
    
    %MS=[3 10 20 9 11 15];% marker size for data 6 points
    MS=[4 8 12 4 8 4];% marker size for data 6 points
    
    %%% 0-1s or 0-3s
    figure(1)
    subplot (3,6, ss+1)
    
   

    %%% First
    tmpdA=PC_A1(x_range2,:);
    tmpdB=PC_B1(x_range2,:);
    tmpdC=PC_C1(x_range2,:);
    tmpdD=PC_D1(x_range2,:);
    
    cTablegreen = [3/255 175/255 122/255];
    cTableOrange = [246/255 170/255 0/255];
    cTableGray = [200/255 200/255 200/255];
    
    if numPC>2
        
        plot3(tmpdA(:,1),tmpdA(:,2),tmpdA(:,3),'-m','Linewidth',1.5);hold on;
        plot3(tmpdB(:,1),tmpdB(:,2),tmpdB(:,3),'-b','Linewidth',1.5);
        plot3(tmpdC(:,1),tmpdC(:,2),tmpdC(:,3),'-','Color',cTablegreen,'Linewidth',1.5);
        plot3(tmpdD(:,1),tmpdD(:,2),tmpdD(:,3),'-','Color',cTableOrange,'Linewidth',1.5);
        
        % t = -1
        plot3(tmpdA(1,1),tmpdA(1,2),tmpdA(1,3),'om','Markerfacecolor','m','Markersize',MS(1));
        plot3(tmpdB(1,1),tmpdB(1,2),tmpdB(1,3),'ob','Markerfacecolor','b','Markersize',MS(1));
        plot3(tmpdC(1,1),tmpdC(1,2),tmpdC(1,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(1));
        plot3(tmpdD(1,1),tmpdD(1,2),tmpdD(1,3),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(1));
        
        % t = 0; odor starts
        plot3(tmpdA(tp1,1),tmpdA(tp1,2),tmpdA(tp1,3),'^m','Markerfacecolor','w','Markersize',MS(2));
        plot3(tmpdB(tp1,1),tmpdB(tp1,2),tmpdB(tp1,3),'^b','Markerfacecolor','w','Markersize',MS(2));
        plot3(tmpdC(tp1,1),tmpdC(tp1,2),tmpdC(tp1,3),'^','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(2));
        plot3(tmpdD(tp1,1),tmpdD(tp1,2),tmpdD(tp1,3),'^','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(2));
        
        % t = 1; odor ends
        plot3(tmpdA(tp2,1),tmpdA(tp2,2),tmpdA(tp2,3),'pm','Markerfacecolor','w','Markersize',MS(3));
        plot3(tmpdB(tp2,1),tmpdB(tp2,2),tmpdB(tp2,3),'pb','Markerfacecolor','w','Markersize',MS(3));
        plot3(tmpdC(tp2,1),tmpdC(tp2,2),tmpdC(tp2,3),'p','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(3));
        plot3(tmpdD(tp2,1),tmpdD(tp2,2),tmpdD(tp2,3),'p','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(3));
        
        % t = 2;
        plot3(tmpdA(tp3,1),tmpdA(tp3,2),tmpdA(tp3,3),'om','Markerfacecolor','m','Markersize',MS(4));
        plot3(tmpdB(tp3,1),tmpdB(tp3,2),tmpdB(tp3,3),'ob','Markerfacecolor','b','Markersize',MS(4));
        plot3(tmpdC(tp3,1),tmpdC(tp3,2),tmpdC(tp3,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(4));
        plot3(tmpdD(tp3,1),tmpdD(tp3,2),tmpdD(tp3,3),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(4));
        
        % t = 3;
        plot3(tmpdA(tp4,1),tmpdA(tp4,2),tmpdA(tp4,3),'sm','Markerfacecolor','w','Markersize',MS(5));
        plot3(tmpdB(tp4,1),tmpdB(tp4,2),tmpdB(tp4,3),'sb','Markerfacecolor','w','Markersize',MS(5));
        plot3(tmpdC(tp4,1),tmpdC(tp4,2),tmpdC(tp4,3),'s','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(5));
        plot3(tmpdD(tp4,1),tmpdD(tp4,2),tmpdD(tp4,3),'s','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(5));
        
        % t = 4;
        plot3(tmpdA(tp5,1),tmpdA(tp5,2),tmpdA(tp5,3),'om','Markerfacecolor','m','Markersize',MS(6));
        plot3(tmpdB(tp5,1),tmpdB(tp5,2),tmpdB(tp5,3),'ob','Markerfacecolor','b','Markersize',MS(6));
        plot3(tmpdC(tp5,1),tmpdC(tp5,2),tmpdC(tp5,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(6));
        plot3(tmpdD(tp5,1),tmpdD(tp5,2),tmpdD(tp5,3),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(6));
        
        xlabel('PC1')
        ylabel('PC2')
        zlabel('PC3')
        axis equal;grid on
%         xlim(xRange)
%         ylim(Index in position 2 is invalid. Array indices must be positive integers or logical values.)
%         zlim([-2.5  2.5])        
    else % when numPC ==2
        

        
        plot(tmpdA(1:tp5,1),tmpdA(1:tp5,2),'-m','Linewidth',2);hold on;
        plot(tmpdB(1:tp5,1),tmpdB(1:tp5,2),'-b','Linewidth',2);
        plot(tmpdC(1:tp5,1),tmpdC(1:tp5,2),'-','Color',cTablegreen,'Linewidth',2);
        plot(tmpdD(1:tp5,1),tmpdD(1:tp5,2),'-','Color',cTableOrange,'Linewidth',2);
        
        % t = -1
        plot(tmpdA(1,1),tmpdA(1,2),'om','Markerfacecolor','m','Markersize',MS(1));
        plot(tmpdB(1,1),tmpdB(1,2),'ob','Markerfacecolor','b','Markersize',MS(1));
        plot(tmpdC(1,1),tmpdC(1,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(1));
        plot(tmpdD(1,1),tmpdD(1,2),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(1));
        
        % t = 0; odor starts
        plot(tmpdA(tp1,1),tmpdA(tp1,2),'^m','Markerfacecolor','w','Markersize',MS(2));
        plot(tmpdB(tp1,1),tmpdB(tp1,2),'^b','Markerfacecolor','w','Markersize',MS(2));
        plot(tmpdC(tp1,1),tmpdC(tp1,2),'^','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(2));
        plot(tmpdD(tp1,1),tmpdD(tp1,2),'^','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(2));
        
        % t = 1; odor ends
        plot(tmpdA(tp2,1),tmpdA(tp2,2),'pm','Markerfacecolor','w','Markersize',MS(3));
        plot(tmpdB(tp2,1),tmpdB(tp2,2),'pb','Markerfacecolor','w','Markersize',MS(3));
        plot(tmpdC(tp2,1),tmpdC(tp2,2),'p','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(3));
        plot(tmpdD(tp2,1),tmpdD(tp2,2),'p','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(3));
        
        % t = 2;
        plot(tmpdA(tp3,1),tmpdA(tp3,2),'om','Markerfacecolor','m','Markersize',MS(4));
        plot(tmpdB(tp3,1),tmpdB(tp3,2),'ob','Markerfacecolor','b','Markersize',MS(4));
        plot(tmpdC(tp3,1),tmpdC(tp3,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(4));
        plot(tmpdD(tp3,1),tmpdD(tp3,2),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(4));
        
        % t = 3;
        plot(tmpdA(tp4,1),tmpdA(tp4,2),'sm','Markerfacecolor','w','Markersize',MS(5));
        plot(tmpdB(tp4,1),tmpdB(tp4,2),'sb','Markerfacecolor','w','Markersize',MS(5));
        plot(tmpdC(tp4,1),tmpdC(tp4,2),'s','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(5));
        plot(tmpdD(tp4,1),tmpdD(tp4,2),'s','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(5));
        
        % t = 4;
        plot(tmpdA(tp5,1),tmpdA(tp5,2),'om','Markerfacecolor','m','Markersize',MS(6));
        plot(tmpdB(tp5,1),tmpdB(tp5,2),'ob','Markerfacecolor','b','Markersize',MS(6));
        plot(tmpdC(tp5,1),tmpdC(tp5,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(6));
        plot(tmpdD(tp5,1),tmpdD(tp5,2),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(6));
        
        xlabel('PC1')
        ylabel('PC2')
        title(titleText{ss})
        axis equal;
        xlim(xRange)
        ylim(yRange)
        set(gcf,'renderer','Painters')
        box off;
        

        
       %% animation
       plotPCAanimation(tmpdA, tmpdB, tmpdC, tmpdD, cTablegreen, cTableOrange, MS, tp1, tp2, tp3, tp4,ss);
        
    end
    
    
    %% color map of distance
    tmp=[ECdist.AB1 ECdist.AC1 ECdist.AD1 ECdist.BC1 ECdist.BD1 ECdist.CD1];
    mMin2=min(min(tmp));
    mMax2=max(max(tmp));
    figure(1)
    subplot(3,6,ss+6+1)

    imagesc(tmp', [0 3]);hold on;
    colormap parula
    set(gca,'ytick',1:6,'yticklabel',{'A-B','A-1','A-2','B-1','B-2','1-2'})
    set(gca,'xtick',[1 21 41 61 81 101],'xticklabel',{'-1','0','1','2','3','4'})
    xlabel('Time from odor onset (s)')
    xlim([1 101])
    plot([21 21],[0 7],'w:')
    plot([41 41],[0 7],'w:')
    plot([81 81],[0 7],'w:')
    set(gcf,'renderer','Painters')
    %colorbar
    
    
    afterOdor500_original.AB =  mean(ECdist.AB1(31:50));
    afterOdor500_original.AC =  mean(ECdist.AC1(31:50));
    afterOdor500_original.AD =  mean(ECdist.AD1(31:50));
    afterOdor500_original.BC =  mean(ECdist.BC1(31:50));
    afterOdor500_original.BD =  mean(ECdist.BD1(31:50));
    afterOdor500_original.CD =  mean(ECdist.CD1(31:50));
 
    
    dist_delay.AB =  mean(ECdist.AB1(51:80));
    dist_delay.AC =  mean(ECdist.AC1(51:80));
    dist_delay.AD =  mean(ECdist.AD1(51:80));
    dist_delay.BC =  mean(ECdist.BC1(51:80));
    dist_delay.BD =  mean(ECdist.BD1(51:80));
    dist_delay.CD =  mean(ECdist.CD1(51:80));
    
    dist_choice.AB =  mean(ECdist.AB1(81:120));
    dist_choice.AC =  mean(ECdist.AC1(81:120));
    dist_choice.AD =  mean(ECdist.AD1(81:120));
    dist_choice.BC =  mean(ECdist.BC1(81:120));
    dist_choice.BD =  mean(ECdist.BD1(81:120));
    dist_choice.CD =  mean(ECdist.CD1(81:120));
    

    
        %% separation index
    ninetyfivePercentile =meanFivePercentile(ss,2);
        
   % SI = (distance - shuffle)/Shuffle; 
   % SI = 1 if identical representation;
   % SI >0 if significant separation
   
    SI.AB (ss) =(afterOdor500_original.AB - ninetyfivePercentile)/ninetyfivePercentile; 
    SI.AC (ss) =(afterOdor500_original.AC - ninetyfivePercentile)/ninetyfivePercentile; 
    SI.AD (ss) =(afterOdor500_original.AD - ninetyfivePercentile)/ninetyfivePercentile; 
    SI.BC (ss) =(afterOdor500_original.BC - ninetyfivePercentile)/ninetyfivePercentile; 
    SI.BD (ss) =(afterOdor500_original.BD - ninetyfivePercentile)/ninetyfivePercentile; 
    SI.CD (ss) =(afterOdor500_original.CD - ninetyfivePercentile)/ninetyfivePercentile; 
    
    
    
    
           %% Integration-Separation index (ISI)
    
   % ISI = (shuffle - distance)/Shuffle; 
   % ISI <0 means significant separation
   % ISI >0 means identical representation;
   
    ISI.AB (ss,1) =(- afterOdor500_original.AB + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.AC (ss,1) =(- afterOdor500_original.AC + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.AD (ss,1) =(- afterOdor500_original.AD + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.BC (ss,1) =(- afterOdor500_original.BC + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.BD (ss,1) =(- afterOdor500_original.BD + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.CD (ss,1) =(- afterOdor500_original.CD + ninetyfivePercentile)/ninetyfivePercentile; 
    
    ISI.AB (ss,2) =(- dist_delay.AB + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.AC (ss,2) =(- dist_delay.AC + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.AD (ss,2) =(- dist_delay.AD + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.BC (ss,2) =(- dist_delay.BC + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.BD (ss,2) =(- dist_delay.BD + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.CD (ss,2) =(- dist_delay.CD + ninetyfivePercentile)/ninetyfivePercentile; 
    
    ISI.AB (ss,3) =(- dist_choice.AB + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.AC (ss,3) =(- dist_choice.AC + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.AD (ss,3) =(- dist_choice.AD + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.BC (ss,3) =(- dist_choice.BC + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.BD (ss,3) =(- dist_choice.BD + ninetyfivePercentile)/ninetyfivePercentile; 
    ISI.CD (ss,3) =(- dist_choice.CD + ninetyfivePercentile)/ninetyfivePercentile; 
    
    
    %% plot distance in the 500ms period after odor offset
    
    figure(1)
    subplot(3,6,ss+12+1)
    
    % distance
    % AB, NaN, AC, AD, NaN, BC, BD, NaN, CD (total: 9)
    plottedData = [afterOdor500_original.AB NaN afterOdor500_original.AC afterOdor500_original.AD NaN afterOdor500_original.BC afterOdor500_original.BD NaN afterOdor500_original.CD];
    
    label={'A-B' '' 'A-1' 'A-2' '' 'B-1' 'B-2' '' '1-2'};
    
    colorTable={[0 0 0] NaN cTablegreen cTableOrange NaN cTablegreen cTableOrange NaN cTableGray};
    colorRectangle(plottedData, colorTable,label)
    hold on
    
    
    plot([0 9.5],[ninetyfivePercentile ninetyfivePercentile],'r')
    
    ylabel('distance during cue')
    set(gcf,'renderer','Painters')
if mode ==1
        ylim([0 2.3])
else
        ylim([0 4.3])
end
    
        

        
end

%% integration-separation index

     fig5=figure(5);
     set( fig5,'position',[300,30,500,500]);
     

     p1=plot([2:6],[ISI.AB(:,1)],'o-','Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'Markerfacecolor',[0 0 0],'Markersize',10,'DisplayName','A-B');hold on;box off
     plot([1],[ISI_ABsession(1)],'o-','Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'Markerfacecolor',[0 0 0],'Markersize',10);hold on;box off
     
     
     p2=plot([2:6],[ISI.AC(:,1)],'o-','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',10,'DisplayName','A-1');hold on;box off
     p3=plot([2:6],[ISI.BD(:,1)],'o-','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',10,'DisplayName','B-2');hold on;box off
     
     plot([0 6.5],[0 0],'r')
     xlim([0.5 6.5])
     ylim([-3 1])
     set(gca,'xtick',[1 2 3 4 5 6],'xticklabel',{'T0','T1','T2','T3','T3','T5'})
     ylabel('ISI')
     legend([p1 p2 p3])
     set(gcf,'renderer','Painters')
     
toc
end

function plotPCAanimation(tmpdA, tmpdB, tmpdC, tmpdD, cTableOrange, cTableCyan, MS, tp1, tp2, tp3, tp4,ss)

t=1:101;
t=(t-1)/20 -1 ;

obj = VideoWriter(['PCAanimationT', num2str(ss),'.avi']);
obj.Quality = 100;
obj.FrameRate = 5;
open(obj);





figure(100+ss)


for n=1:101


        plot(tmpdA(1:n,1),tmpdA(1:n,2),'-m','Linewidth',2);hold on;
        plot(tmpdB(1:n,1),tmpdB(1:n,2),'-b','Linewidth',2);
        plot(tmpdC(1:n,1),tmpdC(1:n,2),'-','Color',cTableOrange,'Linewidth',2);
        plot(tmpdD(1:n,1),tmpdD(1:n,2),'-','Color',cTableCyan,'Linewidth',2);
        
        
        title(num2str(t(n), 'time = %4.3f (sec)'));

         if n>=0 && n<21
            text(3, 2.5, "PRE",'Color','k','FontSize',20)
        end
        
        if n>=21 && n<41
            text(3, 2.5, "CUE",'Color','red','FontSize',20)
        end
        
        if n>=41 && n<81
            text(3, 2.5, "DELAY",'Color','blue','FontSize',20)
        end
        
         if n>=81 
            text(2.7, 2.5, "REWARD",'Color','k','FontSize',20)
        end        
        
        
        
        % t = -1
        plot(tmpdA(1,1),tmpdA(1,2),'om','Markerfacecolor','m','Markersize',MS(1));
        plot(tmpdB(1,1),tmpdB(1,2),'ob','Markerfacecolor','b','Markersize',MS(1));
        plot(tmpdC(1,1),tmpdC(1,2),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(1));
        plot(tmpdD(1,1),tmpdD(1,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(1));
        
        text(2.5,-1.1, "A",'Color','m','FontSize',20)
        text(3,-1.1, "B",'Color','b','FontSize',20)
        text(3.5,-1.1, "1",'Color',cTableOrange,'FontSize',20)
        text(4,-1.1, "2",'Color',cTableCyan,'FontSize',20)
        
        
        if n>= 21
        % t = 0; odor starts
        plot(tmpdA(tp1,1),tmpdA(tp1,2),'^m','Markerfacecolor','w','Markersize',MS(2));
        plot(tmpdB(tp1,1),tmpdB(tp1,2),'^b','Markerfacecolor','w','Markersize',MS(2));
        plot(tmpdC(tp1,1),tmpdC(tp1,2),'^','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(2));
        plot(tmpdD(tp1,1),tmpdD(tp1,2),'^','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(2));
        end
        
        if n>= 41
        % t = 1; odor ends
        plot(tmpdA(tp2,1),tmpdA(tp2,2),'pm','Markerfacecolor','w','Markersize',MS(3));
        plot(tmpdB(tp2,1),tmpdB(tp2,2),'pb','Markerfacecolor','w','Markersize',MS(3));
        plot(tmpdC(tp2,1),tmpdC(tp2,2),'p','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(3));
        plot(tmpdD(tp2,1),tmpdD(tp2,2),'p','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(3));
        

        
        end
        
        if n>= 61
        % t = 2;
        plot(tmpdA(tp3,1),tmpdA(tp3,2),'om','Markerfacecolor','m','Markersize',MS(4));
        plot(tmpdB(tp3,1),tmpdB(tp3,2),'ob','Markerfacecolor','b','Markersize',MS(4));
        plot(tmpdC(tp3,1),tmpdC(tp3,2),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(4));
        plot(tmpdD(tp3,1),tmpdD(tp3,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(4));
        end
        
        if n>= 81
        % t = 3;
        plot(tmpdA(tp4,1),tmpdA(tp4,2),'sm','Markerfacecolor','w','Markersize',MS(5));
        plot(tmpdB(tp4,1),tmpdB(tp4,2),'sb','Markerfacecolor','w','Markersize',MS(5));
        plot(tmpdC(tp4,1),tmpdC(tp4,2),'s','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor','w','Markersize',MS(5));
        plot(tmpdD(tp4,1),tmpdD(tp4,2),'s','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(5));
        end
        
        if n== 101
        % t = 4;
        plot(tmpdA(101),tmpdA(101,2),'om','Markerfacecolor','m','Markersize',MS(6));
        plot(tmpdB(101,1),tmpdB(101,2),'ob','Markerfacecolor','b','Markersize',MS(6));
        plot(tmpdC(101,1),tmpdC(101,2),'o','Color',cTableOrange,'MarkerEdgeColor',cTableOrange,'Markerfacecolor',cTableOrange,'Markersize',MS(6));
        plot(tmpdD(101,1),tmpdD(101,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(6));
        end
        
        xlabel('PC1')
        ylabel('PC2')
        
        axis equal;
        box off;
        xlim([-2  4.5])
        ylim([-1.5  3])
        hold off;
        
        f = getframe(gcf);
        writeVideo(obj, f);
        
        %pause(0.1)
end 

obj.close();

end


function colorRectangle(V, colorTable,label)

% bar graph with error bars
% by Kei Igarashi


% use labels and colorTables as followings:
% label={'AON' 'TT' 'OT' 'APC' 'PPC' 'LEC' 'AMYG'};
% colorTable={'b' [1 0.5 0] 'c' 'r' 'g' 'y' 'm'};



% wright color bar graph for each component of vector V
% with Error bar vector E
% size of array 
%100131 igk
% 


vsize=max(size(V));
for k=1:vsize
    if V(k)==0
        V(k)=0.00000001;
           
    end
    
    if V(k)<0
        rectangle('Position',[k-0.4,V(k),0.8,-V(k)],'FaceColor',colorTable{k})
        hold on
    %line([k,k], [V(k)-E(k),V(k)+E(k)],'LineWidth',2,'Color','k') 
    %line([0,vsize+0.5], [0,0],'LineWidth',1,'Color','k') 
        
    elseif isnan(V(k))
    
    else
    rectangle('Position',[k-0.4,0,0.8,V(k)],'FaceColor',colorTable{k})
    hold on
    %line([k,k], [V(k)-E(k),V(k)+E(k)],'LineWidth',2,'Color','k') 
    
    end
   
end
 set(gca, 'XLim',[0.2 vsize+0.5]);
 set(gca, 'XTick',[1:vsize]);
 set(gca, 'XTickLabel',label,'FontSize',12);
end

