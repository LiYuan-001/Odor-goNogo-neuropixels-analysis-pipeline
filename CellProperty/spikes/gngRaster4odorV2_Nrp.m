%% =====================================
% 4 odor go/nogo task ver.
% 05/08/2018 Tomoaki Nakazono
% Original: equalPlot_maxT.m
% put "Nlx_TS.mat" at the same folder.
% Normal usage: Check_OdorResponse('inFile.txt',1, 0)
%
% Modified 180521 Kei Igarashi
%=======================================

%
% C:\Data\TTList.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% The first line specifies the directory to the first session and also the
% name on the t-file list that will be used for all the sessions listed.
% The t-file list must contain the Mclust t-files, one name for each file.
% If the list contains cells that only occur in some of the sessions, these
% cells will be plotted as having zero firing rate when missing. Placefield
% images will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called placeFieldImages.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is modified from Igarashi lab original code
% I tried not to change it unless necessary, such as path, reading in cells
% Modified to fit neuropixels recording
% Li YUAN, Tohoku, 2026-May

function [ax,posAxisPSTH,raster,PSTH_avg,PSTH_se,h2] = gngRaster4odorV2_Nrp(tsp, odor, onset_ts, premsec,postmsec, blockID)
odorNum = length(unique(odor));

for thisOdor = 1:odorNum
    trialNum_thisodor = sum(odor==thisOdor);
    tiralPlotLim = trialNum_thisodor+3;

    trialInd = odor==thisOdor;
    trialNum = sum(trialInd);
    onsetTs = onset_ts(trialInd);

    trialAfterOptFlag = NaN; % leave the handle for future optotag

    for k = 1:trialNum
        piyo = tsp-(onsetTs(k));
        data_total(k).times = piyo((-premsec < piyo) & (piyo < postmsec))'*1000; % in milisecond
    end
    tmpname = strcat('odor',blockID(thisOdor));
    [ax{thisOdor},posAxisPSTH{thisOdor},tmprate2S,tmprateSE,h2{thisOdor}]=RasterHist4odorV6(data_total, 50, 1,tmpname, premsec*10^3, postmsec*10^3,4,thisOdor,trialAfterOptFlag);

    raster{thisOdor}=data_total;
    PSTH_avg{thisOdor}=tmprate2S;
    PSTH_se{thisOdor}=tmprateSE;

end
end


%%
function [ax,a,b,c,h2]=RasterHist4odorV6(data,bin,line_w, name, premsec, postmsec, lownum ,plotnum,trialAfterOptFlag)
%
% data: 1kHz
% bin : msec

%figure
if ~isempty(data)
    %% raster
    subplot(2,lownum,plotnum)

    hold on
    tmpFr=[];
    for i=1:length(data)
        tS=data(i).times;
        tmpFr=tmpFr+length(tS);

        % for j=1:length(tS)
        %     line([tS(j) tS(j)], [i-1 i], 'LineWidth', line_w);
        %     hold on
        % end

        xpoint = [tS;tS];
        ypos = i;
        ypoint = [ypos+zeros(size(tS))-0.3;ypos+zeros(size(tS))+0.3];
        if ~isempty(tS)
            plot(xpoint,ypoint,'k','LineWidth',line_w)
            hold on
        end
    end
    xlim([-premsec postmsec])
    set(gca,'YDir','reverse')
    ylim([0 length(data)+1])

    line([0 0], [0 length(data)], 'Color','k','LineStyle','--');
    line([1000 1000], [0 length(data)], 'Color','k','LineStyle','--');
    line([3000 3000], [0 length(data)], 'Color','k','LineStyle','--');
    %for i=1:5
    %    line([-20 60], [50*i 50*i]);
    %end

    if ~isnan(trialAfterOptFlag)
        line([-1000 4000], [trialAfterOptFlag-1 trialAfterOptFlag-1], 'Color','r','LineStyle','--');
    end

    if isempty(data)
        tmpy=1;
    else
        tmpy=length(data);
    end
    ylim([0 tmpy]);

    %set(gca,'xtick',[-premsec:1000:postmsec],'xticklabel',[-5:1:5])


    %% firing rate

    trialNum=length(data);

    edges=[-premsec:bin:postmsec]; %Defines time windows for calculating rate
    numSpikes=zeros(trialNum,round(1+(premsec+postmsec)/bin));
    for i=1:trialNum
        numSpikes(i,:)=histc(data(i).times,edges);
    end


    binRate = numSpikes/(bin/1000); %'bin' is in ms; this yields # of spikes/second for each bin

    rate=mean(binRate,1);


    ax=subplot(2,lownum,plotnum+lownum);



    % calculate FR for first N trials and last N trials

    % numInitTrials=5; %
    % numLastTrials=5;

    if trialNum>10
        numInitTrials=10; %
        numLastTrials=10;
    else
        numInitTrials=trialNum; %
        numLastTrials=trialNum;
    end

    %  save('numSpikes.mat')
    binRateInitial = numSpikes(1:numInitTrials,:)/(bin/1000);
    binRateEnd = numSpikes(end-numLastTrials+1:end,:)/(bin/1000);
    binRateAll = numSpikes(1:end,:)/(bin/1000);

    rate2(1,:)=mean(binRateInitial,1);
    rate2(2,:)=mean(binRateEnd,1);
    rate2(3,:)=mean(binRateAll,1);

    rateSE(1,:)=std(binRateInitial)./sqrt(numInitTrials);
    rateSE(2,:)=std(binRateEnd)./sqrt(numLastTrials);
    rateSE(3,:)=std(binRateAll)./sqrt(size(numSpikes,1));


    rate2S=smooth_gaussian2(rate2,100,bin,1);
    rateSE=smooth_gaussian2(rateSE,100,bin,1);

    % bar(edges+bin/2,rate, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor',[0.8 0.8 0.8]);


    %set(gca,'xtick',[-premsec:1000:5000],'xticklabel',[-5:1:5])

    hold on
    posAxis = edges+bin/2;

    colorTable=[1 0 0;  0 0 1];
    colorTable2=[1 0.8 0.8;  0.8 0.8 1];

    % draw SE
    for i=1:2

        X = [posAxis,fliplr(posAxis)];
        Y = [rate2S(i,:)+rateSE(i,:),fliplr(rate2S(i,:)-rateSE(i,:))];
        h=fill(X,Y,'r');

        h.FaceColor=colorTable2(i,:);
        h.EdgeColor = 'none';
        set(h,'facealpha',.5)

        hold on


    end
    maxFR=max(max(rate2S+rateSE));
    xlim([-premsec postmsec])
    if maxFR~=0
        ylim([0 maxFR*1.1])
    end

    for i=1:2
        plot(posAxis,rate2S(i,:),'-', 'Color',colorTable(i,:) ,'LineWidth',2)
    end
    title(name,'FontSize',14);

    a=posAxis;
    b=rate2S;
    c=rateSE;

    %% statistics


    % Statistical test
    testbin = 100;

    edges_t=[-premsec:testbin:postmsec];
    psth_t=zeros(trialNum,round(1+(premsec+postmsec)/testbin));
    for i=1:trialNum
        psth_t(i,:)=histc(data(i).times,edges_t);
    end

    % t-test 5%
    psth_t(isnan(psth_t))=0;

    for n=1:size(psth_t,2)

        %h2(n) = ttest2(psth_t(1:N,n),psth_t(end-N+1:end,n));
        [~,h2(n)] = ranksum(psth_t(1:numInitTrials,n),psth_t(end-numLastTrials+1:end,n));

    end

    %pvalue((pvalue(:)== 1) & (mean(psth_t(1:N,:),2) < mean(psth_t(end-N+1:end,:),2))) = -1
    for i=1:length(h2)
        if h2(i) == 1
            if mean(psth_t(1:numInitTrials,i),1) < mean(psth_t(end-numLastTrials+1:end,i),1)

                line([edges_t(i) edges_t(i+1)], [maxFR+1 maxFR+1],'Color','blue', 'LineWidth', 8);

            else
                line([edges_t(i) edges_t(i+1)], [maxFR+1 maxFR+1],'Color','red', 'LineWidth', 8);

            end
        end
    end
else
    ax=subplot(2,lownum,plotnum+lownum);
    a=NaN;b=NaN;c=NaN;
end


end


function [smoothed2]=smooth_gaussian2(data,sigma,bin_width,dimension)
%sigma bin_width in sec (.004 for 4ms, bin_width=.001 for 1000Hz)

% dimension==1: smooth across columns for data in each rows
% dimension==2: smooth across raws for data in each columns

ratio=10*sigma/bin_width;


if dimension==1

    for i=1:size(data,1)
        data2=[repmat(data(i,1),1,ratio), data(i,:), repmat(data(i,size(data,2)),1,ratio)]; % prolong edges
        edges=[-3*sigma:bin_width:3*sigma];
        kernel=normpdf(edges,0,sigma);
        kernel=kernel*bin_width; %multiply by bin width
        center=ceil(length(edges)/2);
        smoothed=conv(data2,kernel);
        smoothed=smoothed(center:numel(data2)+center-1);
        smoothed2(i,:)=smoothed(ratio+1:numel(data2)-ratio); % cut edges to proper length

    end



elseif  dimension==2

    for i=1:size(data,2)
        data2=[repmat(data(1,i),ratio,1); data(:,i); repmat(data(size(data,1),i),ratio,1)];
        edges=[-3*sigma:bin_width:3*sigma];
        kernel=normpdf(edges,0,sigma);
        kernel=kernel*bin_width; %multiply by bin width
        center=ceil(length(edges)/2);
        smoothed=conv(data2,kernel);
        smoothed=smoothed(center:numel(data2)+center-1);
        smoothed2(:,i)=smoothed(ratio+1:numel(data2)-ratio);

    end


end

end
