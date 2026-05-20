
function [ISI_ABsession] = gngSessionSum2_PCA_AB2_v3(CTS_data2, Time, numPC, xRange, yRange, mode,BehavPerformance,meanFivePercentile_AB)   
tic


sigma = 200; %ms
binSize = 50; %ms

CTS_data2(:,1:121) =   smooth_gaussian2(CTS_data2(:,1:121),  sigma,binSize,1);
CTS_data2(:,122:242) = smooth_gaussian2(CTS_data2(:,122:242),sigma,binSize,1);
% CTS_data2(:,243:363) = smooth_gaussian2(CTS_data2(:,243:363),sigma,binSize,1);
% CTS_data2(:,364:484) = smooth_gaussian2(CTS_data2(:,364:484),sigma,binSize,1);


tp1=21;%0s
tp2=41;%1s
tp3=61;%2s
tp4=81;%3s

 %MS=[3 10 20 9 11 15];% marker size for data 6 points
 %MS=[7 15 20 7 15 7];% marker size for data 6 points
MS=[4 8 12 4 8 4];% marker size for data 6 points
 
 
% Time = matrixRange;
%% PCA Sim1


X=CTS_data2;


%% for shuffled data
% for n=1:size(X,1)
%    m = randperm(4);
%    indexMatrix = [121*(m(1)-1)+1:121*m(1),  121*(m(2)-1)+1:121*m(2),  121*(m(3)-1)+1:121*m(3), 121*(m(4)-1)+1:121*m(4)];
% 
%     X(n,:) = X(n,indexMatrix);
% 
% end


%%
[~,score,~,~, explained] = pca(X');
mMin=min(min(score));
mMax=max(max(score));
PCA.Sim1=score;

AB_ExplainedByPC1_PC2=sum(explained(1:2));
AB_numCells = size(X,1);

%% Trajectories in principal component space
PC_A1=score(Time.odorA1,1:numPC); %numPC =2 or 3
PC_B1=score(Time.odorB1,1:numPC);
% PC_C1=score(Time.odorC1,1:numPC);
% PC_D1=score(Time.odorD1,1:numPC);

%% PCA continuous plot
%%% 0-1s or 0-3s
fig=figure(1);
set( fig,'position',[20,20,1600,900]);
subplot (3,6,1)

tmp_time=Time.odorA1;
%%% First
tmpdA=PC_A1(tmp_time,:);
tmpdB=PC_B1(tmp_time,:);
% tmpdC=PC_C1(tmp_time,:);
% tmpdD=PC_D1(tmp_time,:);

cTablegreen = [3/255 175/255 122/255];
cTableCyan = [77/255 200/255 255/255];
cTableOrange = [246/255 170/255 0/255];
cTableGray = [200/255 200/255 200/255];

if numPC>2
    
    plot3(tmpdA(:,1),tmpdA(:,2),tmpdA(:,3),'-m','Linewidth',1.5);hold on;
    plot3(tmpdB(:,1),tmpdB(:,2),tmpdB(:,3),'-b','Linewidth',1.5);
%     plot3(tmpdC(:,1),tmpdC(:,2),tmpdC(:,3),'-','Color',cTablegreen,'Linewidth',1.5);
%     plot3(tmpdD(:,1),tmpdD(:,2),tmpdD(:,3),'-','Color',cTableCyan,'Linewidth',1.5);
    
    % t = -1
    plot3(tmpdA(1,1),tmpdA(1,2),tmpdA(1,3),'om','Markerfacecolor','m','Markersize',MS(1));
    plot3(tmpdB(1,1),tmpdB(1,2),tmpdB(1,3),'ob','Markerfacecolor','b','Markersize',MS(1));
%     plot3(tmpdC(1,1),tmpdC(1,2),tmpdC(1,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(1));
%     plot3(tmpdD(1,1),tmpdD(1,2),tmpdD(1,3),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(1));
    
    % t = 0; odor starts
    plot3(tmpdA(tp1,1),tmpdA(tp1,2),tmpdA(tp1,3),'sm','Markerfacecolor','m','Markersize',MS(2));
    plot3(tmpdB(tp1,1),tmpdB(tp1,2),tmpdB(tp1,3),'sb','Markerfacecolor','b','Markersize',MS(2));
%     plot3(tmpdC(tp1,1),tmpdC(tp1,2),tmpdC(tp1,3),'s','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(2));
%     plot3(tmpdD(tp1,1),tmpdD(tp1,2),tmpdD(tp1,3),'s','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(2));
%     
    % t = 1; odor ends
    plot3(tmpdA(tp2,1),tmpdA(tp2,2),tmpdA(tp2,3),'pm','Markerfacecolor','m','Markersize',MS(3));
    plot3(tmpdB(tp2,1),tmpdB(tp2,2),tmpdB(tp2,3),'pb','Markerfacecolor','b','Markersize',MS(3));
%     plot3(tmpdC(tp2,1),tmpdC(tp2,2),tmpdC(tp2,3),'p','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(3));
%     plot3(tmpdD(tp2,1),tmpdD(tp2,2),tmpdD(tp2,3),'p','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(3));
%     
    % t = 2;
    plot3(tmpdA(tp3,1),tmpdA(tp3,2),tmpdA(tp3,3),'om','Markerfacecolor','m','Markersize',MS(4));
    plot3(tmpdB(tp3,1),tmpdB(tp3,2),tmpdB(tp3,3),'ob','Markerfacecolor','b','Markersize',MS(4));
%     plot3(tmpdC(tp3,1),tmpdC(tp3,2),tmpdC(tp3,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(4));
%     plot3(tmpdD(tp3,1),tmpdD(tp3,2),tmpdD(tp3,3),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(4));
    
    % t = 3;
    plot3(tmpdA(tp4,1),tmpdA(tp4,2),tmpdA(tp4,3),'om','Markerfacecolor','w','Markersize',MS(5));
    plot3(tmpdB(tp4,1),tmpdB(tp4,2),tmpdB(tp4,3),'ob','Markerfacecolor','w','Markersize',MS(5));
%     plot3(tmpdC(tp4,1),tmpdC(tp4,2),tmpdC(tp4,3),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(5));
%     plot3(tmpdD(tp4,1),tmpdD(tp4,2),tmpdD(tp4,3),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(5));
%     
    % t = 4;
    plot3(tmpdA(end,1),tmpdA(end,2),tmpdA(end,3),'xm','Markerfacecolor','m','Markersize',MS(6));
    plot3(tmpdB(end,1),tmpdB(end,2),tmpdB(end,3),'xb','Markerfacecolor','b','Markersize',MS(6));
%     plot3(tmpdC(end,1),tmpdC(end,2),tmpdC(end,3),'x','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(6));
%     plot3(tmpdD(end,1),tmpdD(end,2),tmpdD(end,3),'x','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(6));
    
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    
    if BehavPerformance ==1
    title(['Correct n =' num2str(AB_numCells) 'cells'])
    else
    title(['Error n =' num2str(AB_numCells) 'cells'])    
    end
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
%     plot(tmpdC(:,1),tmpdC(:,2),'-','Color',cTablegreen,'Linewidth',2);
%     plot(tmpdD(:,1),tmpdD(:,2),'-','Color',cTableCyan,'Linewidth',2);
%     
    % t = -1
    plot(tmpdA(1,1),tmpdA(1,2),'om','Markerfacecolor','m','Markersize',MS(1));
    plot(tmpdB(1,1),tmpdB(1,2),'ob','Markerfacecolor','b','Markersize',MS(1));
%     plot(tmpdC(1,1),tmpdC(1,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(1));
%     plot(tmpdD(1,1),tmpdD(1,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(1));
%     
    % t = 0; odor starts
    plot(tmpdA(tp1,1),tmpdA(tp1,2),'^m','Markerfacecolor','w','Markersize',MS(2));
    plot(tmpdB(tp1,1),tmpdB(tp1,2),'^b','Markerfacecolor','w','Markersize',MS(2));
%     plot(tmpdC(tp1,1),tmpdC(tp1,2),'^','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(2));
%     plot(tmpdD(tp1,1),tmpdD(tp1,2),'^','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(2));
    
    % t = 1; odor ends
    plot(tmpdA(tp2,1),tmpdA(tp2,2),'pm','Markerfacecolor','w','Markersize',MS(3));
    plot(tmpdB(tp2,1),tmpdB(tp2,2),'pb','Markerfacecolor','w','Markersize',MS(3));
%     plot(tmpdC(tp2,1),tmpdC(tp2,2),'p','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(3));
%     plot(tmpdD(tp2,1),tmpdD(tp2,2),'p','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(3));
    
    % t = 2;
    plot(tmpdA(tp3,1),tmpdA(tp3,2),'om','Markerfacecolor','m','Markersize',MS(4));
    plot(tmpdB(tp3,1),tmpdB(tp3,2),'ob','Markerfacecolor','b','Markersize',MS(4));
%     plot(tmpdC(tp3,1),tmpdC(tp3,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(4));
%     plot(tmpdD(tp3,1),tmpdD(tp3,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(4));
    
    % t = 3;
    plot(tmpdA(tp4,1),tmpdA(tp4,2),'sm','Markerfacecolor','w','Markersize',MS(5));
    plot(tmpdB(tp4,1),tmpdB(tp4,2),'sb','Markerfacecolor','w','Markersize',MS(5));
%     plot(tmpdC(tp4,1),tmpdC(tp4,2),'s','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor','w','Markersize',MS(5));
%     plot(tmpdD(tp4,1),tmpdD(tp4,2),'s','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor','w','Markersize',MS(5));
    
    % t = 4;
    plot(tmpdA(end,1),tmpdA(end,2),'om','Markerfacecolor','m','Markersize',MS(6));
    plot(tmpdB(end,1),tmpdB(end,2),'ob','Markerfacecolor','b','Markersize',MS(6));
%     plot(tmpdC(end,1),tmpdC(end,2),'o','Color',cTablegreen,'MarkerEdgeColor',cTablegreen,'Markerfacecolor',cTablegreen,'Markersize',MS(6));
%     plot(tmpdD(end,1),tmpdD(end,2),'o','Color',cTableCyan,'MarkerEdgeColor',cTableCyan,'Markerfacecolor',cTableCyan,'Markersize',MS(6));
%     
    xlabel('PC1')
    ylabel('PC2')
    if BehavPerformance ==1
    title(['Correct n =' num2str(AB_numCells) 'cells'])
    else
    title(['Error n =' num2str(AB_numCells) 'cells'])    
    end
    axis equal;
    box off;
    xlim(xRange)
    ylim(yRange)
    
end





  %% plotting Euclidian distances
% 
% fig3=figure(3);
% set( fig3,'position',[200,20,800,300]);




PC_A1=score(Time.odorA1,1:numPC);
PC_B1=score(Time.odorB1,1:numPC);



% % -----------------------------------------------------------------------------
% %%% Euclidean distance

normalization = 0;
PC_C1 = NaN;
PC_D1 = NaN;
odorNumber = 2;
% --------------------------
ECdist = gngSessionSum2_PCA_ECdist_v5(PC_A1, PC_B1, PC_C1, PC_D1, normalization, odorNumber);

% --------------------------



fivePercentile =meanFivePercentile_AB(1);
ninetyfivePercentile =meanFivePercentile_AB(2);
if  mode ==0 || mode ==6
 ninetyfivePercentile = 0.81;
else
  ninetyfivePercentile = 1.2;  
end

tmp=[ECdist.AB1];

    figure(1)
    subplot(3,6,7)
    

imagesc(tmp', [0 3]);hold on;
colormap parula
set(gca,'ytick',1,'yticklabel',{'A-B'})
set(gca,'xtick',[1 21 41 61 81 101],'xticklabel',{'-1','0','1','2','3','4'})
xlabel('Time from odor onset (s)')
  xlim([1 101])
plot([21 21],[0 7],'w:')
plot([41 41],[0 7],'w:')
plot([81 81],[0 7],'w:')
colorbar




%% calculate shuffling data

% mean of the original data

afterOdor500_original.AB =  mean(ECdist.AB1(31:50));
dist_delay.AB =  mean(ECdist.AB1(51:80));
dist_choice.AB =  mean(ECdist.AB1(81:100));


% plot distance in the 500ms period after odor offset




    figure(1)
    subplot(3,6,13)

% distance
% AB, NaN, AC, AD, NaN, BC, BD, NaN, CD (total: 9)
% plottedData = [afterOdor500_original.AB NaN NaN NaN NaN NaN NaN NaN NaN];
% 
% label={'AB' '' '' '' '' '' '' '' ''};
% 
%     colorTable={[0 0 0] NaN cTablegreen cTableCyan NaN cTablegreen cTableCyan NaN cTableGray};
%     colorRectangle(plottedData, colorTable,label)
%     hold on
% 
% 
% plot([0 9.5],[ninetyfivePercentile ninetyfivePercentile],'r')
% 
% ylabel('distance')
% 
% ylim([0 4.3])

% 
% 
% fig=figure(2);
% set( fig,'position',[20,20,1600,900]);
%      subplot(3,6,1)
% 
% % distance
% % AB, NaN, AC, AD, NaN, BC, BD, NaN, CD (total: 9)
% plottedData = [dist_delay.AB NaN NaN NaN NaN NaN NaN NaN NaN];
% 
% label={'AB' '' '' '' '' '' '' '' ''};
%     colorTable={[0 0 0] NaN cTablegreen cTableCyan NaN cTablegreen cTableCyan NaN cTableGray};
%     colorRectangle(plottedData, colorTable,label)
%     hold on
% 
% plot([0 9.5],[ninetyfivePercentile ninetyfivePercentile],'r')
% ylabel('delay distance')
% ylim([0 4.3])
% 
% 

     % subplot(3,6,7)

% distance
% AB, NaN, AC, AD, NaN, BC, BD, NaN, CD (total: 9)
plottedData = [dist_choice.AB NaN NaN NaN NaN NaN NaN NaN NaN];

label={'A-B' '' '' '' '' '' '' '' ''};
    colorTable={[0 0 0] NaN cTablegreen cTableOrange NaN cTablegreen cTableOrange NaN cTableGray};
    colorRectangle(plottedData, colorTable,label)
    hold on

plot([0 9.5],[ninetyfivePercentile ninetyfivePercentile],'r')
ylabel('choice distance')
ylim([0 4.3])

ISI_ABsession(1) =(- afterOdor500_original.AB + ninetyfivePercentile)/ninetyfivePercentile; 


ISI_ABsession(2) =(- dist_delay.AB + ninetyfivePercentile)/ninetyfivePercentile; 
ISI_ABsession(3) =(- dist_choice.AB + ninetyfivePercentile)/ninetyfivePercentile; 


toc
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
