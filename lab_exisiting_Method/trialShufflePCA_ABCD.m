function [meanFivePercentile] = trialShufflePCA_ABCD(folderName,odorIndex_allCells, rst_all_allCells, Time, x_range, numPC, mode, BehavPerformance,inhStimTag,firstLast,shuffle)

%calculate 5 percentile of shuffled distance -500 to +500 ms around odor ending
% 200820 Kei

bin=0.05;%sec, bin for psth
bin_size=50;%ms, bin for psth
win=[-2 6];%sec, bin for raster and psth(8 sec, 161 bins data)
showFig =0;

%% make distination for saving folder
    cd(folderName);
    
    mkdir ABCD_shuffle
    folder =strcat(folderName, '\ABCD_shuffle');
    cd (folder)
    

%% 1. generating 1000 shuffled data

if shuffle == 0

for m= 1:1000 % shuffle 1000 times
    
    for ss= 1:5 % T1 - T5
        % calculate PSTH for each cell
        for i= 1:size(odorIndex_allCells,2)
            
            odorIndex = odorIndex_allCells{1,i};
            
            %assing shuffled trials
            odorIndex_shuffle =  odorIndex(randperm(length(odorIndex)));
            
            rst_all = rst_all_allCells{1,i};
            for s= 1:4 % A, B, C, D
                trialID=find(odorIndex_shuffle==s);
                numTrials=size(trialID,1);
                increment =round((numTrials-10)/4);
                
                rst=[];
                
                if inhStimTag==0 && mode ==1
                    trialRange(1,:) = numTrials-9:numTrials;
                    trialRange(2,:) = numTrials-12:numTrials-3;
                    trialRange(3,:) = round(numTrials/2)-4:round(numTrials/2)+5;
                    trialRange(4,:) = increment:increment+9;
                    trialRange(5,:) = 1:10;
                    
                    
                else
                    
                    trialRange(1,:) = 1:10; % first10
                    trialRange(2,:) = increment:increment+9; %T2
                    %trialRange(3,:) = increment*2:increment*2+9; %T3
                    trialRange(3,:) = round(numTrials/2)-4:round(numTrials/2)+5;  % middle10
                    
                    %trialRange(4,:) = numTrials-10:numTrials-1;   %T4: for Sim1 good
                    trialRange(4,:) = numTrials-12:numTrials-3;   %T4: for Sim1 bad
                    %trialRange(4,:) = increment*3:increment*3+9;  %T4
                    trialRange(5,:) = numTrials-9:numTrials; % last10
                    
                    
                end
                
                
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
        
        CTS_data2(:,1:121) =   smooth_gaussian2(CTS_data2(:,1:121),  sigma,binSize,1);
        CTS_data2(:,122:242) = smooth_gaussian2(CTS_data2(:,122:242),sigma,binSize,1);
        CTS_data2(:,243:363) = smooth_gaussian2(CTS_data2(:,243:363),sigma,binSize,1);
        CTS_data2(:,364:484) = smooth_gaussian2(CTS_data2(:,364:484),sigma,binSize,1);
        
        
        
        %% PCA
        [~,score,~,~, ~] = pca(CTS_data2');
        
        
        PC_A1=score(Time.odorA1,1:numPC); %numPC =2 or 3
        PC_B1=score(Time.odorB1,1:numPC);
        PC_C1=score(Time.odorC1,1:numPC);
        PC_D1=score(Time.odorD1,1:numPC);
        
        % ExplainedByPC1_numPC=sum(explained(1:numPC))
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
        
%         odor1000.AB(ss) =  mean(ECdist.AB1(31:50));
%         odor1000.AC(ss) =  mean(ECdist.AC1(31:50));
%         odor1000.AD(ss) =  mean(ECdist.AD1(31:50));
%         odor1000.BC(ss) =  mean(ECdist.BC1(31:50));
%         odor1000.BD(ss) =  mean(ECdist.BD1(31:50));
%         odor1000.CD(ss) =  mean(ECdist.CD1(31:50));
        
         % mean of 1- 5 s after odor onset
        odor1000.AB(ss) =  mean(ECdist.AB1(41:120));
        odor1000.AC(ss) =  mean(ECdist.AC1(41:120));
        odor1000.AD(ss) =  mean(ECdist.AD1(41:120));
        odor1000.BC(ss) =  mean(ECdist.BC1(41:120));
        odor1000.BD(ss) =  mean(ECdist.BD1(41:120));
        odor1000.CD(ss) =  mean(ECdist.CD1(41:120));
        
        
        %% PCA continuous plot
        
        if showFig ==1
            
            mMin=min(min(score));
            mMax=max(max(score));
            
            tp1=21;%0s
            tp2=41;%1s
            tp3=61;%2s
            tp4=81;%3s
            
            %MS=[3 10 20 9 11 15];% marker size for data 6 points
            MS=[7 15 20 7 15 7];% marker size for data 6 points
            
            %%% 0-1s or 0-3s
            close all
            fig=figure(1);
            set( fig,'position',[200,20,800,600]);
            
            tmp_time=Time.time1;
            %%% First
            tmpdA=PC_A1(tmp_time,:);
            tmpdB=PC_B1(tmp_time,:);
            tmpdC=PC_C1(tmp_time,:);
            tmpdD=PC_D1(tmp_time,:);
            
            cTablegreen = [3/255 175/255 122/255];
            cTableCyan = [77/255 200/255 255/255];
            cTableGray = [200/255 200/255 200/255];
            
            if numPC>2
                
                plot3(tmpdA(:,1),tmpdA(:,2),tmpdA(:,3),'-m','Linewidth',1.5);hold on;
                plot3(tmpdB(:,1),tmpdB(:,2),tmpdB(:,3),'-b','Linewidth',1.5);
                plot3(tmpdC(:,1),tmpdC(:,2),tmpdC(:,3),'-','Color',cTablegreen,'Linewidth',1.5);
                plot3(tmpdD(:,1),tmpdD(:,2),tmpdD(:,3),'-','Color',cTableCyan,'Linewidth',1.5);
                
                % t = -1
                plot3(tmpdA(1,1),tmpdA(1,2),tmpdA(1,3),'om','Markerfacecolor','m','Markersize',MS(1));
                plot3(tmpdB(1,1),tmpdB(1,2),tmpdB(1,3),'ob','Markerfacecolor','b','Markersize',MS(1));
                plot3(tmpdC(1,1),tmpdC(1,2),tmpdC(1,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(1));
                plot3(tmpdD(1,1),tmpdD(1,2),tmpdD(1,3),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(1));
                
                % t = 0; odor starts
                plot3(tmpdA(tp1,1),tmpdA(tp1,2),tmpdA(tp1,3),'sm','Markerfacecolor','m','Markersize',MS(2));
                plot3(tmpdB(tp1,1),tmpdB(tp1,2),tmpdB(tp1,3),'sb','Markerfacecolor','b','Markersize',MS(2));
                plot3(tmpdC(tp1,1),tmpdC(tp1,2),tmpdC(tp1,3),'s','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(2));
                plot3(tmpdD(tp1,1),tmpdD(tp1,2),tmpdD(tp1,3),'s','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(2));
                
                % t = 1; odor ends
                plot3(tmpdA(tp2,1),tmpdA(tp2,2),tmpdA(tp2,3),'pm','Markerfacecolor','m','Markersize',MS(3));
                plot3(tmpdB(tp2,1),tmpdB(tp2,2),tmpdB(tp2,3),'pb','Markerfacecolor','b','Markersize',MS(3));
                plot3(tmpdC(tp2,1),tmpdC(tp2,2),tmpdC(tp2,3),'p','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(3));
                plot3(tmpdD(tp2,1),tmpdD(tp2,2),tmpdD(tp2,3),'p','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(3));
                
                % t = 2;
                plot3(tmpdA(tp3,1),tmpdA(tp3,2),tmpdA(tp3,3),'om','Markerfacecolor','m','Markersize',MS(4));
                plot3(tmpdB(tp3,1),tmpdB(tp3,2),tmpdB(tp3,3),'ob','Markerfacecolor','b','Markersize',MS(4));
                plot3(tmpdC(tp3,1),tmpdC(tp3,2),tmpdC(tp3,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(4));
                plot3(tmpdD(tp3,1),tmpdD(tp3,2),tmpdD(tp3,3),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(4));
                
                % t = 3;
                plot3(tmpdA(tp4,1),tmpdA(tp4,2),tmpdA(tp4,3),'om','Markerfacecolor','w','Markersize',MS(5));
                plot3(tmpdB(tp4,1),tmpdB(tp4,2),tmpdB(tp4,3),'ob','Markerfacecolor','w','Markersize',MS(5));
                plot3(tmpdC(tp4,1),tmpdC(tp4,2),tmpdC(tp4,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(5));
                plot3(tmpdD(tp4,1),tmpdD(tp4,2),tmpdD(tp4,3),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(5));
                
                % t = 4;
                plot3(tmpdA(end,1),tmpdA(end,2),tmpdA(end,3),'xm','Markerfacecolor','m','Markersize',MS(6));
                plot3(tmpdB(end,1),tmpdB(end,2),tmpdB(end,3),'xb','Markerfacecolor','b','Markersize',MS(6));
                plot3(tmpdC(end,1),tmpdC(end,2),tmpdC(end,3),'x','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(6));
                plot3(tmpdD(end,1),tmpdD(end,2),tmpdD(end,3),'x','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(6));
                
                xlabel('PC1')
                ylabel('PC2')
                zlabel('PC3')
                axis equal;grid on
                
            else % when numPC ==2
                
                %     plot(tmpdA(1:tp2-10,1),tmpdA(1:tp2-10,2),'-m','Linewidth',1);hold on;
                %     plot(tmpdB(1:tp2-10,1),tmpdB(1:tp2-10,2),'-b','Linewidth',1);
                %     plot(tmpdC(1:tp2-10,1),tmpdC(1:tp2-10,2),'-','Color',cTablegreen,'Linewidth',1);
                %     plot(tmpdD(1:tp2-10,1),tmpdD(1:tp2-10,2),'-','Color',cTableCyan,'Linewidth',1);
                %
                %     plot(tmpdA(tp2-10:tp2+10,1),tmpdA(tp2-10:tp2+10,2),'-m','Linewidth',2);hold on;
                %     plot(tmpdB(tp2-10:tp2+10,1),tmpdB(tp2-10:tp2+10,2),'-b','Linewidth',2);
                %     plot(tmpdC(tp2-10:tp2+10,1),tmpdC(tp2-10:tp2+10,2),'-','Color',cTablegreen,'Linewidth',2);
                %     plot(tmpdD(tp2-10:tp2+10,1),tmpdD(tp2-10:tp2+10,2),'-','Color',cTableCyan,'Linewidth',2);
                %
                %     plot(tmpdA(tp2+10:end,1),tmpdA(tp2+10:end,2),'-m','Linewidth',1);hold on;
                %     plot(tmpdB(tp2+10:end,1),tmpdB(tp2+10:end,2),'-b','Linewidth',1);
                %     plot(tmpdC(tp2+10:end,1),tmpdC(tp2+10:end,2),'-','Color',cTablegreen,'Linewidth',1);
                %     plot(tmpdD(tp2+10:end,1),tmpdD(tp2+10:end,2),'-','Color',cTableCyan,'Linewidth',1);
                
                plot(tmpdA(:,1),tmpdA(:,2),'-m','Linewidth',2);hold on;
                plot(tmpdB(:,1),tmpdB(:,2),'-b','Linewidth',2);
                plot(tmpdC(:,1),tmpdC(:,2),'-','Color',cTablegreen,'Linewidth',2);
                plot(tmpdD(:,1),tmpdD(:,2),'-','Color',cTableCyan,'Linewidth',2);
                
                % t = -1
                plot(tmpdA(1,1),tmpdA(1,2),'om','Markerfacecolor','m','Markersize',MS(1));
                plot(tmpdB(1,1),tmpdB(1,2),'ob','Markerfacecolor','b','Markersize',MS(1));
                plot(tmpdC(1,1),tmpdC(1,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(1));
                plot(tmpdD(1,1),tmpdD(1,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(1));
                
                % t = 0; odor starts
                plot(tmpdA(tp1,1),tmpdA(tp1,2),'^m','Markerfacecolor','w','Markersize',MS(2));
                plot(tmpdB(tp1,1),tmpdB(tp1,2),'^b','Markerfacecolor','w','Markersize',MS(2));
                plot(tmpdC(tp1,1),tmpdC(tp1,2),'^','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(2));
                plot(tmpdD(tp1,1),tmpdD(tp1,2),'^','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(2));
                
                % t = 1; odor ends
                plot(tmpdA(tp2,1),tmpdA(tp2,2),'pm','Markerfacecolor','w','Markersize',MS(3));
                plot(tmpdB(tp2,1),tmpdB(tp2,2),'pb','Markerfacecolor','w','Markersize',MS(3));
                plot(tmpdC(tp2,1),tmpdC(tp2,2),'p','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(3));
                plot(tmpdD(tp2,1),tmpdD(tp2,2),'p','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(3));
                
                % t = 2;
                plot(tmpdA(tp3,1),tmpdA(tp3,2),'om','Markerfacecolor','m','Markersize',MS(4));
                plot(tmpdB(tp3,1),tmpdB(tp3,2),'ob','Markerfacecolor','b','Markersize',MS(4));
                plot(tmpdC(tp3,1),tmpdC(tp3,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(4));
                plot(tmpdD(tp3,1),tmpdD(tp3,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(4));
                
                % t = 3;
                plot(tmpdA(tp4,1),tmpdA(tp4,2),'sm','Markerfacecolor','w','Markersize',MS(5));
                plot(tmpdB(tp4,1),tmpdB(tp4,2),'sb','Markerfacecolor','w','Markersize',MS(5));
                plot(tmpdC(tp4,1),tmpdC(tp4,2),'s','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(5));
                plot(tmpdD(tp4,1),tmpdD(tp4,2),'s','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(5));
                
                % t = 4;
                plot(tmpdA(end,1),tmpdA(end,2),'om','Markerfacecolor','m','Markersize',MS(6));
                plot(tmpdB(end,1),tmpdB(end,2),'ob','Markerfacecolor','b','Markersize',MS(6));
                plot(tmpdC(end,1),tmpdC(end,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(6));
                plot(tmpdD(end,1),tmpdD(end,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(6));
                
                xlabel('PC1')
                ylabel('PC2')
                
                axis equal;
                box off;
                xlim([-2  4.5])
                ylim([-1.5  3])
                
            end
 
        end

    end
    
           outf=['ABCD_PCAdist_shuffle_',num2str(m)];
           save(outf,'odor1000') 
    % if mode==0
    %     if  BehavPerformance==1
    %         outf=['N:\soma\populations\4odorsTaskResults\summaryJasonSim1\shuffle\PCAdist_shuffle_',num2str(m)];
    %         save(outf,'odor1000')
    %     else
    %         outf=['N:\soma\populations\4odorsTaskResults\summaryJasonSim1\shuffle_bad\PCAdist_shuffle_',num2str(m)];
    %         save(outf,'odor1000')
    %     end
    % elseif mode==1
    %     if inhStimTag==1 %inhibited
    %         outf=['N:\soma\populations\4odorsTaskResults\spike_data\pend_Sim1xDAT4\shuffle\PCAdist_shuffle_',num2str(m)];
    %         save(outf,'odor1000')
    %     else
    %         outf=['N:\soma\populations\4odorsTaskResults\spike_data\pend_Sim1xDAT4\shuffle_nonInhibited\PCAdist_shuffle_',num2str(m)];
    %         save(outf,'odor1000')
    %     end
    % elseif mode==4 % DAT inh Correct
    %         outf=['N:\soma\populations\4odorsTaskResults\spike_data\pend_Sim1xDAT4\shuffleCorrect\PCAdist_shuffle_',num2str(m)];
    %         save(outf,'odor1000')                
    % 
    % end
end

elseif shuffle == 1
    disp('Shuffle is skipped. The shuffle value is calculated using the previous shuffle result.')

end


%% 2. pull out 95th percentile

cd(folder);

for m= 1:1000
   
    outf=['ABCD_PCAdist_shuffle_',num2str(m)];
load(outf,'odor1000')

% for m= 1:1000
%     if mode==0 || mode ==6 || mode ==7|| mode ==8|| mode ==9|| mode ==2|| mode ==3
%         if  BehavPerformance==1
%             outf=['N:\soma\populations\4odorsTaskResults\summaryJasonSim1\shuffle\PCAdist_shuffle_',num2str(m)];
%         else
%             outf=['N:\soma\populations\4odorsTaskResults\summaryJasonSim1\shuffle_bad\PCAdist_shuffle_',num2str(m)];
%         end
%     elseif mode==1|| mode ==5 
%         if inhStimTag==1 %inhibited
%             outf=['N:\soma\populations\4odorsTaskResults\spike_data\pend_Sim1xDAT4\shuffle\PCAdist_shuffle_',num2str(m)];
%         else
%            outf=['N:\soma\populations\4odorsTaskResults\spike_data\pend_Sim1xDAT4\shuffle_nonInhibited\PCAdist_shuffle_',num2str(m)];
%         end
%      elseif mode==4 % DAT inh Correct
%             outf=['N:\soma\populations\4odorsTaskResults\spike_data\pend_Sim1xDAT4\shuffleCorrect\PCAdist_shuffle_',num2str(m)];
% 
%     end
    
    
% load(outf)
        
        
dist.AB(:,m)=odor1000.AB;
dist.AC(:,m)=odor1000.AC;
dist.AD(:,m)=odor1000.AD;
dist.BC(:,m)=odor1000.BC;
dist.BD(:,m)=odor1000.BD;
dist.CD(:,m)=odor1000.CD;
end



for ss = 1:5 %T1 - T5
    
    sortedAB = sort(dist.AB(ss,:));
    sortedAC = sort(dist.AC(ss,:));
    sortedAD = sort(dist.AD(ss,:));
    sortedBC = sort(dist.BC(ss,:));
    sortedBD = sort(dist.BD(ss,:));
    sortedCD = sort(dist.CD(ss,:));
    
    fivePercentile(:,1) = [sortedAB(50) sortedAB(950)];
    fivePercentile(:,2) = [sortedAC(50) sortedAC(950)];
    fivePercentile(:,3) = [sortedAD(50) sortedAD(950)];
    fivePercentile(:,4) = [sortedBC(50) sortedBC(950)];
    fivePercentile(:,5) = [sortedBD(50) sortedBD(950)];
    fivePercentile(:,6) = [sortedCD(50) sortedCD(950)];
    
    
    if mode ==1 || inhStimTag ==0
        
        meanFivePercentile(ss,:) = [min(fivePercentile(1,:)) min(fivePercentile(2,:))];
        
    else
        meanFivePercentile(ss,:) = mean(fivePercentile,2);
        

    end
    
    
    
end
%% plotting shuffled data - use good last10

% figure
% hist(sortedBC, 30);
% hold on
% box off
% plot([meanFivePercentile(2) meanFivePercentile(2)], [0 120],'r')
% 


end

